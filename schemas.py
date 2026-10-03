from typing import Literal, Optional

from pydantic import BaseModel


class Normalized(BaseModel):
    gene: Optional[str] = None
    transcript: Optional[str] = None
    hgvs_c: Optional[str] = None
    hgvs_p: Optional[str] = None
    rsid: Optional[str] = None
    vcf: Optional[str] = None


class Evidence(BaseModel):
    source: Literal["VEP", "ClinVar", "gnomAD"]
    finding: str
    url: Optional[str] = None


class Conflict(BaseModel):
    type: str
    detail: str


class Verdict(BaseModel):
    variant_input: str
    normalized: Normalized
    verdict: Literal[
        "Pathogenic",
        "Likely pathogenic",
        "Uncertain significance",
        "Likely benign",
        "Benign",
        "Insufficient evidence",
    ]
    confidence: Literal["High", "Medium", "Low"]
    summary: str
    evidence: list[Evidence]
    conflicts: list[Conflict]
    caveats: list[str]
