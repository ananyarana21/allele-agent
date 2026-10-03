import re
import threading
import time
import urllib.parse

import requests

ENSEMBL = "https://rest.ensembl.org"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GNOMAD = "https://gnomad.broadinstitute.org/api"
ANCESTRY_GROUPS = ("afr", "amr", "ami", "asj", "eas", "fin", "mid", "nfe", "sas")
SINGLE_BASE_CHROM = re.compile(r"^(\d+|X|Y|MT)-(\d+)-([ACGTN]+)-([ACGTN]+)$")

memo = {}
ncbi_lock = threading.Lock()
ncbi_last = [0.0]


def http(method, url, **kwargs):
    error = "request failed"
    for attempt in range(2):
        try:
            if "eutils.ncbi" in url:
                with ncbi_lock:
                    wait = 0.4 - (time.time() - ncbi_last[0])
                    if wait > 0:
                        time.sleep(wait)
                    ncbi_last[0] = time.time()
            response = requests.request(method, url, timeout=15, **kwargs)
            if response.status_code == 429 or response.status_code >= 500:
                error = f"HTTP {response.status_code}"
                time.sleep(1.5)
                continue
            try:
                body = response.json()
            except ValueError:
                return None, f"HTTP {response.status_code}: response was not JSON"
            if response.status_code >= 400:
                message = body.get("error") if isinstance(body, dict) else None
                return None, message or f"HTTP {response.status_code}"
            return body, None
        except requests.RequestException as exc:
            error = type(exc).__name__
            time.sleep(1)
    return None, f"{error} after one retry"


def memoized(key, fn):
    if key in memo:
        return memo[key]
    result = fn()
    if "error" not in result:
        memo[key] = result
    return result


def trim(d):
    return {k: v for k, v in d.items() if v not in (None, "", [], {})}


def reference_window(chrom, start, end):
    body, error = http(
        "GET",
        f"{ENSEMBL}/sequence/region/human/{chrom}:{start}..{end}:1",
        params={"content-type": "application/json"},
    )
    if error:
        return None
    return body["seq"].upper()


def left_align(chrom, pos, ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return pos, ref, alt
    start = max(1, pos - 80)
    seq = reference_window(chrom, start, pos)
    if not seq:
        return pos, ref, alt
    for _ in range(80):
        if ref and alt and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
        elif (not ref or not alt) and pos - 1 >= start:
            pos -= 1
            base = seq[pos - start]
            ref, alt = base + ref, base + alt
        else:
            break
    return pos, ref, alt


def spdi_to_vcf(spdi):
    accession, pos0, deleted, inserted = spdi.split(":")
    match = re.match(r"NC_0+(\d+)\.\d+", accession)
    if not match:
        return None
    chrom = {"23": "X", "24": "Y"}.get(match.group(1), match.group(1))
    pos, ref, alt = int(pos0) + 1, deleted, inserted
    if not ref or not alt or len(ref) != 1 or len(alt) != 1:
        pos, ref, alt = left_align(chrom, pos, ref, alt)
    return f"{chrom}-{pos}-{ref}-{alt}"


def ensembl_frequencies(rsid):
    body, error = http(
        "GET",
        f"{ENSEMBL}/variation/human/{rsid}",
        params={"content-type": "application/json", "pops": "1"},
    )
    if error:
        return None, error
    rows = [
        p for p in body.get("populations", [])
        if p["population"].startswith(("gnomADe:", "gnomADg:"))
    ]
    return rows, None


def pick_allele(candidates, rsid):
    if len(candidates) == 1 or not rsid:
        return candidates[0]
    rows, error = ensembl_frequencies(rsid)
    if error or not rows:
        return candidates[0]

    def score(candidate):
        alt = candidate["vcf"].split("-")[3]
        values = [r["frequency"] for r in rows if r["allele"] == alt]
        overall = [r["frequency"] for r in rows if r["allele"] == alt and r["population"].endswith(":ALL")]
        return max(overall) if overall else (max(values) if values else 0)

    return max(candidates, key=score)


def normalize_raw(variant):
    text = variant.strip()
    url = f"{ENSEMBL}/variant_recoder/human/{urllib.parse.quote(text, safe=':>_.*()=')}"
    body, error = http(
        "GET", url,
        params={"content-type": "application/json", "fields": "hgvsg,id", "vcf_string": "1"},
    )
    if error:
        return {"error": f"Ensembl Variant Recoder failed: {error}. Ask the user for a different format (rsID, HGVS with transcript, or genomic HGVS).", "source_url": url}
    candidates = []
    for entry in body if isinstance(body, list) else []:
        for allele, info in entry.items():
            if not isinstance(info, dict):
                continue
            hgvsg = next((h for h in info.get("hgvsg", []) if h.startswith("NC_")), None)
            vcf = next((v for v in info.get("vcf_string", []) if SINGLE_BASE_CHROM.match(v)), None)
            if hgvsg and vcf:
                rsids = [i for i in info.get("id", []) if re.match(r"^rs\d+$", i)]
                candidates.append({"allele": allele, "hgvs_g": hgvsg, "vcf": vcf, "rsids": rsids})
    if not candidates:
        return {"error": "Ensembl Variant Recoder returned no genomic mapping for this input. Ask the user for a different format.", "source_url": url}
    rsid = text.lower() if re.match(r"^rs\d+$", text.lower()) else (candidates[0]["rsids"] or [None])[0]
    chosen = pick_allele(candidates, rsid)
    chrom, pos, ref, alt = SINGLE_BASE_CHROM.match(chosen["vcf"]).groups()
    pos, ref, alt = left_align(chrom, int(pos), ref, alt)
    result = {
        "input": text,
        "rsid": rsid,
        "hgvs_g": chosen["hgvs_g"],
        "vcf": f"{chrom}-{pos}-{ref}-{alt}",
        "other_alleles": [c["vcf"] for c in candidates if c is not chosen] if len(candidates) > 1 else None,
        "allele_choice_note": "rsID has several alt alleles; the one with the highest gnomAD frequency was analysed" if len(candidates) > 1 else None,
        "source_url": url,
    }
    annotation = vep_raw(result["hgvs_g"])
    if "error" not in annotation:
        result.update(trim({
            "gene": annotation.get("gene"),
            "transcript": annotation.get("transcript"),
            "hgvs_c": annotation.get("hgvs_c"),
            "hgvs_p": annotation.get("hgvs_p"),
        }))
    return trim(result)


def normalize(variant):
    return memoized(("normalize", variant.strip()), lambda: normalize_raw(variant))


def vep_raw(hgvs_g):
    return memoized(("vep", hgvs_g), lambda: vep_fetch(hgvs_g))


def vep_fetch(hgvs_g):
    url = f"{ENSEMBL}/vep/human/hgvs/{urllib.parse.quote(hgvs_g, safe=':>_.')}"
    body, error = http(
        "GET", url,
        params={"content-type": "application/json", "SIFT": "b", "PolyPhen": "b", "canonical": "1",
                "hgvs": "1", "numbers": "1", "CADD": "1", "REVEL": "1", "mane": "1"},
    )
    if error:
        return {"error": f"Ensembl VEP failed: {error}", "source_url": url}
    if not isinstance(body, list) or not body:
        return {"error": "Ensembl VEP returned no annotation", "source_url": url}
    entry = body[0]
    consequences = entry.get("transcript_consequences", [])
    chosen = (
        next((t for t in consequences if t.get("mane_select")), None)
        or next((t for t in consequences if t.get("canonical") and t.get("gene_symbol") and t.get("biotype") == "protein_coding"), None)
        or next((t for t in consequences if t.get("gene_symbol")), None)
    )
    result = {"most_severe_consequence": entry.get("most_severe_consequence"), "source_url": url}
    if not chosen:
        result["note"] = "No gene transcript consequence (variant may be intergenic)"
        return trim(result)
    hgvs_c = chosen.get("hgvsc")
    result.update({
        "gene": chosen.get("gene_symbol"),
        "transcript": hgvs_c.split(":")[0] if hgvs_c else chosen.get("transcript_id"),
        "mane_select": chosen.get("mane_select"),
        "impact": chosen.get("impact"),
        "transcript_consequences": chosen.get("consequence_terms"),
        "hgvs_c": hgvs_c,
        "hgvs_p": chosen.get("hgvsp"),
        "amino_acids": chosen.get("amino_acids"),
        "protein_position": chosen.get("protein_start"),
        "exon": chosen.get("exon"),
        "intron": chosen.get("intron"),
        "sift": trim({"term": chosen.get("sift_prediction"), "score": chosen.get("sift_score")}),
        "polyphen": trim({"term": chosen.get("polyphen_prediction"), "score": chosen.get("polyphen_score")}),
        "cadd_phred": chosen.get("cadd_phred"),
        "revel": chosen.get("revel"),
    })
    return trim(result)


def normalize_variant(variant):
    return normalize(variant)


def get_vep_annotation(variant):
    normalized = normalize(variant)
    if "error" in normalized:
        return normalized
    return vep_raw(normalized["hgvs_g"])


def clinvar_candidates(term):
    url = f"{EUTILS}/esearch.fcgi"
    found, error = http("GET", url, params={"db": "clinvar", "term": term, "retmode": "json", "retmax": 25})
    if error:
        return None, 0, error
    result = found["esearchresult"]
    ids = result.get("idlist", [])
    if not ids:
        return {}, int(result.get("count", 0)), None
    summary, error = http(
        "GET", f"{EUTILS}/esummary.fcgi",
        params={"db": "clinvar", "id": ",".join(ids), "retmode": "json"},
    )
    if error:
        return None, 0, error
    records = [summary["result"][i] for i in ids if i in summary.get("result", {})]
    return records, int(result["count"]), None


def get_clinvar_raw(variant):
    normalized = normalize(variant)
    if "error" in normalized:
        return normalized
    rsid = normalized.get("rsid")
    terms = ([f"{rsid}[RS]"] if rsid else []) + [normalized["hgvs_g"]]
    records, total = [], 0
    for term in terms:
        records, total, error = clinvar_candidates(term)
        if error:
            return {"error": f"ClinVar (NCBI E-utilities) failed: {error}"}
        if records:
            break
    if not records:
        return {"found": False, "note": "No ClinVar record found for this variant", "source_url": f"https://www.ncbi.nlm.nih.gov/clinvar/?term={urllib.parse.quote(terms[0])}"}
    singles = [r for r in records if len(r.get("variation_set", [])) == 1]
    pool = singles or records
    rs_number = rsid[2:] if rsid else None
    match = "first_hit_unconfirmed"
    chosen = pool[0]
    by_rsid = [r for r in pool if any(x.get("db_source") == "dbSNP" and x.get("db_id") == rs_number for v in r["variation_set"] for x in v.get("variation_xrefs", []))]
    exact = None
    for record in (by_rsid or pool)[:6]:
        spdi = record["variation_set"][0].get("canonical_spdi")
        if spdi and spdi_to_vcf(spdi) == normalized["vcf"]:
            exact = record
            break
    if exact:
        chosen, match = exact, "same_allele"
    elif by_rsid:
        chosen, match = by_rsid[0], "same_rsid_allele_unconfirmed"
    germline = chosen.get("germline_classification") or {}
    legacy = chosen.get("clinical_significance") or {}
    significance = germline.get("description") or legacy.get("description")
    review = germline.get("review_status") or legacy.get("review_status")
    last = germline.get("last_evaluated") or legacy.get("last_evaluated")
    traits = [t.get("trait_name") for t in germline.get("trait_set", []) if t.get("trait_name")]
    variation_id = chosen["uid"]
    return trim({
        "found": True,
        "variation_id": variation_id,
        "title": chosen.get("title"),
        "clinical_significance": significance,
        "review_status": review,
        "num_submitters": len(chosen.get("supporting_submissions", {}).get("scv", [])) or None,
        "conditions": traits[:5],
        "last_evaluated": last.split(" ")[0] if last else None,
        "record_match": match,
        "total_hits": total,
        "source_url": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{variation_id}/",
    })


def get_clinvar(variant):
    return memoized(("clinvar", variant.strip()), lambda: get_clinvar_raw(variant))


def gnomad_query(variant_id):
    query = """query($v:String!){variant(variantId:$v,dataset:gnomad_r4){
        exome{ac an af homozygote_count populations{id ac an}}
        genome{ac an af homozygote_count populations{id ac an}}}}"""
    body, error = http("POST", GNOMAD, json={"query": query, "variables": {"v": variant_id}})
    if error:
        return None, error
    variant = (body.get("data") or {}).get("variant")
    if variant:
        return variant, None
    messages = [e.get("message", "") for e in body.get("errors", [])]
    if all("not found" in m.lower() for m in messages):
        return {}, None
    return None, "; ".join(messages) or "empty response"


def summarize_gnomad(variant, variant_id):
    ac = an = hom = 0
    groups = {}
    for key in ("exome", "genome"):
        part = variant.get(key)
        if not part:
            continue
        ac += part.get("ac") or 0
        an += part.get("an") or 0
        hom += part.get("homozygote_count") or 0
        for pop in part.get("populations", []):
            if pop["id"] in ANCESTRY_GROUPS:
                g = groups.setdefault(pop["id"], [0, 0])
                g[0] += pop["ac"] or 0
                g[1] += pop["an"] or 0
    rates = {k: v[0] / v[1] for k, v in groups.items() if v[1]}
    top = max(rates, key=rates.get) if rates else None
    return trim({
        "found": True,
        "dataset": "gnomAD v4 (exomes + genomes)",
        "allele_count": ac,
        "allele_number": an,
        "global_af": float(f"{ac / an:.4g}") if an else None,
        "homozygote_count": hom,
        "max_population": top,
        "max_population_af": float(f"{rates[top]:.4g}") if top else None,
        "source_url": f"https://gnomad.broadinstitute.org/variant/{variant_id}?dataset=gnomad_r4",
    })


def ensembl_fallback(normalized, reason):
    rsid = normalized.get("rsid")
    if not rsid:
        return None
    rows, error = ensembl_frequencies(rsid)
    if error:
        return None
    alt = normalized["vcf"].split("-")[3]
    rows = [r for r in rows if r["allele"] == alt]
    if not rows:
        return None
    top = max(rows, key=lambda r: r["frequency"])
    overall = next((r for r in rows if r["population"] == "gnomADe:ALL"), None)
    return trim({
        "found": True,
        "dataset": "gnomAD via Ensembl (fallback)",
        "fallback_reason": reason,
        "global_af": float(f"{overall['frequency']:.4g}") if overall else None,
        "max_population": top["population"],
        "max_population_af": float(f"{top['frequency']:.4g}"),
        "source_url": f"https://rest.ensembl.org/variation/human/{rsid}?pops=1;content-type=application/json",
    })


def get_population_frequency_raw(variant):
    normalized = normalize(variant)
    if "error" in normalized:
        return normalized
    variant_id = normalized["vcf"]
    gnomad, error = gnomad_query(variant_id)
    if gnomad:
        return summarize_gnomad(gnomad, variant_id)
    fallback = ensembl_fallback(normalized, error or "not found in gnomAD")
    if fallback:
        return fallback
    if error:
        return {"error": f"gnomAD failed: {error}"}
    return {
        "found": False,
        "note": "Not observed in gnomAD v4 (very rare or never seen)",
        "queried_variant_id": variant_id,
        "source_url": f"https://gnomad.broadinstitute.org/variant/{variant_id}?dataset=gnomad_r4",
    }


def get_population_frequency(variant):
    return memoized(("gnomad", variant.strip()), lambda: get_population_frequency_raw(variant))


FUNCTIONS = {
    "normalize_variant": normalize_variant,
    "get_vep_annotation": get_vep_annotation,
    "get_clinvar": get_clinvar,
    "get_population_frequency": get_population_frequency,
}

DESCRIPTIONS = {
    "normalize_variant": "Call this first. Resolves the user's variant (rsID or HGVS) to rsID, genomic HGVS, VCF-style ID, gene, transcript, HGVS c. and p.",
    "get_vep_annotation": "Get Ensembl VEP functional consequence, impact, SIFT, PolyPhen, CADD and REVEL for the variant. Call after normalize_variant.",
    "get_clinvar": "Get the ClinVar clinical significance, review status, submitter count and conditions for the variant. Call after normalize_variant.",
    "get_population_frequency": "Get gnomAD population allele frequencies, popmax and homozygote counts for the variant. Call after normalize_variant.",
}

TOOL_SCHEMAS = [
    {
        "name": name,
        "description": description,
        "input_schema": {
            "type": "object",
            "properties": {"variant": {"type": "string", "description": "The variant exactly as the user entered it"}},
            "required": ["variant"],
        },
    }
    for name, description in DESCRIPTIONS.items()
]


def run_tool(name, variant):
    function = FUNCTIONS.get(name)
    if not function:
        return {"error": f"Unknown tool {name}"}
    try:
        return function(variant)
    except Exception as exc:
        return {"error": f"{name} crashed: {type(exc).__name__}: {exc}"}
