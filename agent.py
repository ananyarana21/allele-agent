import json
import os
import time

from dotenv import load_dotenv
from google import genai
from google.genai import errors, types
from pydantic import ValidationError

from schemas import Verdict
from tools import TOOL_SCHEMAS, run_tool

load_dotenv()

MODEL = os.environ.get("GEMINI_MODEL", "gemini-3.8-flash")
MAX_ITERATIONS = 8
TRUNCATE_CHARS = 1500

SYSTEM_PROMPT = """You are a variant interpretation assistant for educational use. You are not a clinician.

Process:
1. Call normalize_variant first.
2. Then call get_vep_annotation, get_clinvar and get_population_frequency using the normalized forms. Call them in parallel when possible.
3. Merge the evidence. Never state a fact that is not in a tool result. If a tool returned nothing or errored, say so.
4. Detect conflicts explicitly. Flag when:
   - ClinVar is Pathogenic or Likely pathogenic but SIFT and PolyPhen both predict tolerated/benign
   - ClinVar is Benign or Likely benign but predictors say deleterious/damaging
   - Allele frequency above 1% in any gnomAD population while ClinVar or predictors suggest pathogenic
   - ClinVar has conflicting interpretations or low review status (no assertion criteria, single submitter)
   - Variant absent from ClinVar and gnomAD (insufficient evidence, not proof of pathogenicity)
5. Give a verdict: one of Pathogenic, Likely pathogenic, Uncertain significance, Likely benign, Benign, Insufficient evidence. Base it on the evidence, and when ClinVar has a well-reviewed classification, lean on it and say so.
6. Explain in plain English at a high-school reading level. Short sentences.
7. Finish with one JSON object matching the schema, nothing after it.

Always include the disclaimer that this is not medical advice."""

SCHEMA_INSTRUCTION = (
    "Final answer format: one JSON object, nothing else, with exactly these keys: "
    "variant_input (string), "
    "normalized (object with gene, transcript, hgvs_c, hgvs_p, rsid, vcf; each string or null), "
    "verdict (one of Pathogenic, Likely pathogenic, Uncertain significance, Likely benign, Benign, Insufficient evidence), "
    "confidence (High, Medium or Low), "
    "summary (string, plain English, ends with a not-medical-advice disclaimer), "
    "evidence (list of objects with source = VEP, ClinVar or gnomAD, finding, url or null), "
    "conflicts (list of objects with type and detail; empty list if none), "
    "caveats (list of strings)."
)


def get_client():
    key = os.environ.get("GEMINI_API_KEY") or os.environ.get("GOOGLE_API_KEY")
    if not key:
        try:
            import streamlit as st
            key = st.secrets.get("GEMINI_API_KEY")
        except Exception:
            key = None
    if not key:
        raise RuntimeError("GEMINI_API_KEY is not set (environment variable or Streamlit secret)")
    return genai.Client(api_key=key)


def extract_json(text):
    start, end = text.find("{"), text.rfind("}")
    if start == -1 or end <= start:
        raise ValueError("no JSON object found")
    return json.loads(text[start:end + 1])


def parse_verdict(text):
    return Verdict.model_validate(extract_json(text))


def generate(client, contents, config):
    for attempt in range(4):
        try:
            return client.models.generate_content(model=MODEL, contents=contents, config=config)
        except errors.APIError as exc:
            if attempt == 3 or exc.code not in (429, 500, 503):
                raise
            time.sleep(5 * 2 ** attempt)


def run_agent(variant, on_event, execute=run_tool, enabled_tools=None):
    client = get_client()
    declarations = [
        types.FunctionDeclaration(name=t["name"], description=t["description"], parameters_json_schema=t["input_schema"])
        for t in TOOL_SCHEMAS
        if enabled_tools is None or t["name"] in enabled_tools
    ]
    system = SYSTEM_PROMPT + "\n\n" + SCHEMA_INSTRUCTION
    if enabled_tools is not None:
        missing = [t["name"] for t in TOOL_SCHEMAS if t["name"] not in enabled_tools]
        system += f"\n\nThese tools are unavailable in this run: {', '.join(missing)}. Do not call them, and note their absence in caveats."
    base = dict(
        system_instruction=system,
        tools=[types.Tool(function_declarations=declarations)],
        automatic_function_calling=types.AutomaticFunctionCallingConfig(disable=True),
    )
    config = types.GenerateContentConfig(**base)
    no_tools = types.GenerateContentConfig(
        **base,
        tool_config=types.ToolConfig(function_calling_config=types.FunctionCallingConfig(mode="NONE")),
    )
    contents = [types.Content(role="user", parts=[types.Part.from_text(text=f"Interpret this variant: {variant}")])]
    final_text = ""
    for _ in range(MAX_ITERATIONS):
        response = generate(client, contents, config)
        calls = response.function_calls or []
        if not calls:
            final_text = response.text or ""
            contents.append(response.candidates[0].content)
            break
        contents.append(response.candidates[0].content)
        parts = []
        for call in calls:
            arguments = dict(call.args or {})
            on_event({"type": "tool_call", "name": call.name, "input": arguments})
            started = time.time()
            output = execute(call.name, arguments.get("variant", variant))
            latency = int((time.time() - started) * 1000)
            serialized = json.dumps(output)
            on_event({
                "type": "tool_result",
                "name": call.name,
                "result": output,
                "truncated": serialized[:TRUNCATE_CHARS] + ("..." if len(serialized) > TRUNCATE_CHARS else ""),
                "latency_ms": latency,
            })
            parts.append(types.Part.from_function_response(name=call.name, response=output))
        contents.append(types.Content(role="user", parts=parts))
    else:
        contents.append(types.Content(role="user", parts=[types.Part.from_text(text="Stop calling tools. Give the final JSON answer now.")]))
        response = generate(client, contents, no_tools)
        final_text = response.text or ""
        contents.append(response.candidates[0].content)
    try:
        verdict = parse_verdict(final_text)
    except (ValueError, ValidationError) as exc:
        contents.append(types.Content(role="user", parts=[types.Part.from_text(text=f"That was not valid. Error: {exc}. Reply with valid JSON only, matching the schema.")]))
        final_text = generate(client, contents, no_tools).text or ""
        verdict = parse_verdict(final_text)
    on_event({"type": "final", "verdict": verdict, "raw": final_text})
    return verdict
