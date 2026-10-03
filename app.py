import json

import pandas as pd
import streamlit as st

from agent import run_agent
from tools import run_tool

EXAMPLES = {
    "rs80357906 (BRCA1)": "rs80357906",
    "rs334 (HBB sickle cell)": "rs334",
    "rs113993960 (CFTR F508del)": "rs113993960",
    "rs1801133 (MTHFR C677T)": "rs1801133",
}

BADGE_COLORS = {
    "Pathogenic": ("#c0392b", "#ffffff"),
    "Likely pathogenic": ("#e67e22", "#ffffff"),
    "Uncertain significance": ("#7f8c8d", "#ffffff"),
    "Likely benign": ("#a9dfbf", "#1b4332"),
    "Benign": ("#27ae60", "#ffffff"),
    "Insufficient evidence": ("#bdc3c7", "#2c3e50"),
}

st.set_page_config(page_title="Allele Agent", page_icon="🧬", layout="wide")


@st.cache_data(ttl=3600, show_spinner=False)
def cached_tool(name, variant):
    result = run_tool(name, variant)
    if "error" in result:
        raise RuntimeError(json.dumps(result))
    return result


def execute(name, variant):
    try:
        return cached_tool(name, variant)
    except RuntimeError as exc:
        return json.loads(str(exc))


def set_variant(value):
    st.session_state["variant"] = value


def render_verdict(verdict, raw):
    background, foreground = BADGE_COLORS[verdict.verdict]
    st.markdown(
        f"<div style='display:inline-block;padding:0.6rem 1.4rem;border-radius:0.6rem;"
        f"background:{background};color:{foreground};font-size:2rem;font-weight:700'>"
        f"{verdict.verdict}</div>"
        f"<span style='margin-left:1rem;font-size:1.1rem'>Confidence: <b>{verdict.confidence}</b></span>",
        unsafe_allow_html=True,
    )
    normalized = {k: v for k, v in verdict.normalized.model_dump().items() if v}
    if normalized:
        st.caption("  ·  ".join(f"{k}: {v}" for k, v in normalized.items()))
    st.write(verdict.summary)
    if verdict.conflicts:
        st.warning("**Conflicts and flags**\n\n" + "\n".join(f"- **{c.type}**: {c.detail}" for c in verdict.conflicts))
    if verdict.evidence:
        st.subheader("Evidence")
        table = pd.DataFrame([{"Source": e.source, "Finding": e.finding, "Link": e.url} for e in verdict.evidence])
        st.dataframe(
            table,
            hide_index=True,
            use_container_width=True,
            column_config={"Link": st.column_config.LinkColumn("Link", display_text="open")},
        )
    if verdict.caveats:
        st.subheader("Caveats")
        for caveat in verdict.caveats:
            st.markdown(f"- {caveat}")
    with st.expander("Raw agent output"):
        st.code(raw, language="json")


def analyze(variant):
    final = {}
    with st.status("Agent working...", expanded=True) as status:

        def on_event(event):
            if event["type"] == "tool_call":
                st.write(f"🔧 **{event['name']}**  `{json.dumps(event['input'])}`")
            elif event["type"] == "tool_result":
                flag = "❌" if "error" in event["result"] else "✅"
                st.write(f"{flag} **{event['name']}** returned in {event['latency_ms']} ms")
                with st.expander(f"{event['name']} result", expanded=False):
                    st.json(event["result"])
            else:
                final.update(event)

        try:
            run_agent(variant, on_event, execute=execute)
        except Exception as exc:
            status.update(label="Agent failed", state="error")
            st.error(f"The agent failed: {type(exc).__name__}: {exc}")
            return None
        status.update(label="Agent finished", state="complete", expanded=False)
    return final


def main():
    st.title("🧬 Allele Agent")
    st.caption("An LLM agent that calls Ensembl, ClinVar and gnomAD, merges the evidence and flags conflicts.")
    st.info("Educational use only. This is not medical advice and not a clinical-grade interpretation.")
    left, right = st.columns([3, 1])
    with left:
        st.text_input(
            "Variant",
            key="variant",
            placeholder="rs80357906, NM_007294.4:c.5266dup, ENST00000357654.9:c.5266dup or 17:g.43057062dup",
        )
    with right:
        st.write("")
        st.write("")
        run = st.button("Analyze", type="primary", use_container_width=True)
    columns = st.columns(len(EXAMPLES))
    for column, (label, value) in zip(columns, EXAMPLES.items()):
        column.button(label, on_click=set_variant, args=(value,), use_container_width=True)
    if not run:
        return
    variant = st.session_state.get("variant", "").strip()
    if not variant:
        st.error("Please enter a variant first.")
        return
    with st.spinner("Resolving variant..."):
        probe = execute("normalize_variant", variant)
    if "error" in probe:
        st.error(f"Could not resolve this variant. {probe['error']}")
        return
    final = analyze(variant)
    if final:
        render_verdict(final["verdict"], final["raw"])


main()
