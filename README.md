# Allele Agent

Allele Agent explains the effect of a human genetic variant. Paste a variant and an LLM agent decides which public genomics APIs to call, combines the evidence from Ensembl VEP, ClinVar and gnomAD, flags conflicts between sources, and returns a plain-English verdict with links. The Streamlit UI shows the live tool-call trace so you can see what the agent did.

The tool loop is written by hand on the Google Gemini SDK (`google-genai`). No LangChain or LangGraph.

## Files

- `app.py`: Streamlit UI
- `agent.py`: agent loop and JSON repair
- `tools.py`: API tools (Ensembl Variant Recoder, VEP, ClinVar, gnomAD)
- `schemas.py`: `Verdict` schema

## Accepted inputs

- rsID: `rs80357906`
- HGVS with transcript: `NM_007294.4:c.5266dup`
- HGVS genomic: `NC_000017.11:g.43057065dup`

## Setup

```
pip install -r requirements.txt
```

Set a Gemini API key (free key at https://aistudio.google.com/apikey), either in a `.env` file or as an environment variable:

```
GEMINI_API_KEY=your-key
```

Optionally set `GEMINI_MODEL` to use a different model. Do not commit your key.

## Run

```
streamlit run app.py
```

## Notes

- Multi-allelic rsIDs: the allele with the highest gnomAD frequency is analysed and the others are listed.
- Indels are left-aligned before the gnomAD lookup.
- API calls have a 15 second timeout and one retry. Failures are passed to the agent as errors rather than raised.
