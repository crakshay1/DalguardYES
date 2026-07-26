from __future__ import annotations

import json
import random
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, List, Tuple

import pandas as pd
import plotly.express as px
import streamlit as st

import orbs_duplex as orbs
import candidates as candgen
import riboguard_ga_engine_clean as ga


# -----------------------------------------------------------------------------
# Page setup and styling
# -----------------------------------------------------------------------------

st.set_page_config(
    page_title="Shine DalguardYES",
    page_icon="",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown(
    """
    <style>
    header[data-testid="stHeader"] {
        display: none;
    }

    .stAppToolbar {
        display: none;
    }

    footer {
        display: none;
    }

    #MainMenu {
        visibility: hidden;
    }

    [data-testid="stSidebar"] {
        transform: none !important;
        margin-left: 0 !important;
        width: 320px !important;
        min-width: 320px !important;
        max-width: 320px !important;
        visibility: visible !important;
        opacity: 1 !important;
    }

    button[aria-label="Collapse sidebar"],
    button[title="Collapse sidebar"],
    [data-testid="collapsedControl"],
    [data-testid="stSidebarCollapseButton"] {
        display: none !important;
    }

    * {
        border-radius: 0 !important;
    }

    .stApp {
        background: #0e1117;
        color: #f5f5f5;
    }

    [data-testid="stSidebar"] {
        background: #0b0c10;
        border-right: 1px solid #b30000;
        color: #f5f5f5;
    }

    .block-container {
        padding-top: 1.8rem;
        padding-bottom: 2rem;
        background: #0e1117;
    }

    .hero-title {
        font-size: 2.1rem;
        font-weight: 800;
        line-height: 1.05;
        margin-bottom: 0.15rem;
        color: #ff3333;
    }

    .hero-subtitle {
        color: #cccccc;
        margin-bottom: 1rem;
    }

    .metric-card {
        padding: 1rem 1.05rem;
        border: 1px solid #b30000;
        border-radius: 0px;
        background: #111316;
        box-shadow: none;
        min-height: 92px;
    }

    .metric-label {
        color: #ff9999;
        font-size: 0.86rem;
        margin-bottom: 0.35rem;
    }

    .metric-value {
        color: #ff3333;
        font-size: 1.65rem;
        font-weight: 800;
        margin-bottom: 0;
    }

    .metric-value-blue,
    .metric-value-purple {
        color: #ff3333;
    }

    .panel-title {
        color: #ff3333;
        font-weight: 750;
        font-size: 1.05rem;
        margin: 0.25rem 0 0.55rem 0;
    }

    .small-note {
        color: #999999;
        font-size: 0.82rem;
    }

    .pill-good {
        display: inline-block;
        color: #ffffff;
        background: #b30000;
        border-radius: 0px;
        padding: 0.18rem 0.55rem;
        font-weight: 800;
        font-size: 0.78rem;
    }

    .sequence-box {
        font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", monospace;
        padding: 0.75rem;
        border: 1px solid #b30000;
        border-radius: 0px;
        background: #111316;
        overflow-x: auto;
        white-space: nowrap;
        font-size: 0.92rem;
        color: #f5f5f5;
    }

    .rbs-span,
    .aug-span {
        background: #b30000;
        color: #ffffff;
        padding: 0.12rem 0.22rem;
        border-radius: 0px;
        font-weight: 800;
    }

    div[data-testid="stMetric"] {
        background: #111316;
        border: 1px solid #b30000;
        padding: 0.6rem 0.75rem;
        border-radius: 0px;
        color: #f5f5f5;
    }

    .stButton>button,
    button,
    .stButton button {
        border-radius: 0 !important;
        background: #b30000 !important;
        color: #ffffff !important;
        border: 1px solid #b30000 !important;
    }
    </style>
    """,
    unsafe_allow_html=True,
)


# -----------------------------------------------------------------------------
# Pipeline helpers
# -----------------------------------------------------------------------------


def make_seed_with_orbs_duplex(
    *,
    search_sequence: str,
    five_prime_flank: str,
    cds_start: str,
    name: str,
    asd_tail_nt: int,
    window_k: int,
    orth_asd: str,
    temp_c: float,
) -> Tuple[Dict[str, Any], pd.DataFrame, pd.DataFrame]:
    """Run the ORBS-duplex seed step in memory, without intermediate files."""
    full_seq = orbs.normalize_rna(search_sequence)
    mrna5 = orbs.normalize_rna(five_prime_flank)
    mrna3 = orbs.normalize_rna(cds_start)

    if len(full_seq) < window_k:
        raise ValueError("ORBS search sequence must be at least as long as window_k.")

    tail_nt = max(1, min(int(asd_tail_nt), len(full_seq)))
    start_1based = len(full_seq) - tail_nt + 1
    end_1based = len(full_seq)
    seq_window = full_seq[start_1based - 1:end_1based]

    hits = orbs.window_scan(seq_window, start_1based, int(window_k), float(temp_c))
    merged = orbs.merge_top_hits(hits, seq_window, start_1based, float(temp_c))

    if not hits or not merged:
        raise RuntimeError("ORBS-duplex did not find a usable seed region.")

    best = merged[0]
    core = hits[0]

    seed = {
        "name": orbs.sanitize_name(name),
        "five_prime_flank": mrna5,
        "canonical_rbs": "",
        "sequence":orth_asd,
        "rbs": orbs.revcomp_rna(best.subseq),
        "core": orbs.revcomp_rna(core.subseq),
        "spacer": "",
        "cds_start": mrna3,
        "mutable_regions": ["rbs", "five_prime_flank", "spacer"],
    }

    hits_df = pd.DataFrame([
        {
            "abs_start": h.abs_start,
            "abs_end": h.abs_end,
            "subseq": h.subseq,
            "dg_self_duplex": h.dg,
            "seed_core_rbs": orbs.revcomp_rna(h.subseq),
        }
        for h in hits
    ])

    merged_df = pd.DataFrame([
        {
            "abs_start": m.abs_start,
            "abs_end": m.abs_end,
            "length": len(m.subseq),
            "subseq": m.subseq,
            "dg_self_duplex": m.dg,
            "seed_rbs": orbs.revcomp_rna(m.subseq),
        }
        for m in merged
    ])

    return seed, hits_df, merged_df


@st.cache_data(show_spinner=False)
def run_full_pipeline_cached(
    search_sequence: str,
    orth_asd: str,
    wt_asd: str,
    five_prime_flank: str,
    cds_start: str,
    name: str,
    asd_tail_nt: int,
    window_k: int,
    temp_c: float,
    candidates_per_seed: int,
    standby_start: int,
    spacer_min: int,
    spacer_max: int,
    max_tries: int,
    generations: int,
    population_size: int,
    elite_fraction: float,
    wt_penalty_constant: float,
    rng_seed: int,
    top_n: int,
) -> Dict[str, Any]:
    """ORBS-duplex -> candidates.py -> GA optimizer -> dashboard dataset."""
    random.seed(int(rng_seed))

    seed_record, hits_df, merged_df = make_seed_with_orbs_duplex(
        search_sequence=search_sequence,
        five_prime_flank=five_prime_flank,
        cds_start=cds_start,
        name=name,
        asd_tail_nt=asd_tail_nt,
        window_k=window_k,
        orth_asd=orth_asd,
        temp_c=temp_c,
    )

    initial_candidates = candgen.generate_candidates(
        seed=seed_record,
        n=int(candidates_per_seed),
        standby_start=int(standby_start),
        spacer_len_min=int(spacer_min),
        spacer_len_max=int(spacer_max),
        max_tries=int(max_tries),
    )

    if not initial_candidates:
        raise RuntimeError(
            "candidates.py produced zero candidates. Try increasing max tries or changing spacer settings."
        )

    for c in initial_candidates:
        if "id" not in c:
            c["id"] = c.get("name", f"seed_{len(c)}")

    evals, history, binding_sites, final_population = ga.run_ga(
        initial_candidates=initial_candidates,
        orth_anti_sd=orth_asd,
        wt_anti_sd=wt_asd,
        default_flank=five_prime_flank,
        default_cds_start=cds_start,
        generations=int(generations),
        population_size=int(population_size),
        elite_fraction=float(elite_fraction),
        wt_penalty_constant=float(wt_penalty_constant),
        seed=int(rng_seed),
    )

    dataset = ga.build_dashboard_dataset(
        evals=evals,
        history=history,
        binding_sites=binding_sites,
        inputs={
            "orthogonalAntiSD": orth_asd,
            "wtAntiSD": wt_asd,
            "cdsStart": cds_start,
            "targetExpression": "High",
            "wtPenaltyConstant": wt_penalty_constant,
        },
        top_n=int(top_n),
    )

    eval_rows = [asdict(e) for e in evals]
    binding_rows = [asdict(s) for s in binding_sites]

    dataset["pipeline"] = {
        "orbsSeed": seed_record,
        "orbsHits": hits_df.to_dict(orient="records"),
        "orbsMergedHits": merged_df.to_dict(orient="records"),
        "initialCandidates": initial_candidates,
        "finalPopulation": final_population,
        "steps": [
            "ORBS-duplex seed selection",
            "candidates.py initial candidate generation",
            "GA optimization with ΔG/TIR/T-score objective",
        ],
    }

    dataset["allCandidateRows"] = eval_rows
    dataset["allBindingSiteRows"] = binding_rows

    return dataset


def to_df(records: List[Dict[str, Any]]) -> pd.DataFrame:
    return pd.DataFrame(records) if records else pd.DataFrame()


def fmt_sci(x: Any) -> str:
    try:
        x = float(x)
        if abs(x) >= 1e4 or (abs(x) < 1e-2 and x != 0):
            return f"{x:.2e}"
        return f"{x:.4g}"
    except Exception:
        return "—"


def highlight_sequence(
    full_seq: str,
    structure: str,
    rbs_start: int,
    rbs_end: int,
    aug_start: int,
    aug_end: int,
) -> str:
    before = full_seq[:rbs_start]
    rbs = full_seq[rbs_start:rbs_end]
    middle = full_seq[rbs_end:aug_start]
    aug = full_seq[aug_start:aug_end]
    after = full_seq[aug_end:]

    before_struct = structure[:rbs_start]
    rbs_struct = structure[rbs_start:rbs_end]
    middle_struct = structure[rbs_end:aug_start]
    aug_struct = structure[aug_start:aug_end]
    after_struct = structure[aug_end:]

    return (
        f"<div class='sequence-box'><strong>Sequence:</strong> 5′-"
        f"{before}<span class='rbs-span'>{rbs}</span>{middle}"
        f"<span class='aug-span'>{aug}</span>{after}-3′</div>"
        f"<div class='sequence-box'><strong>Structure:</strong> "
        f"{before_struct}<span class='rbs-span'>{rbs_struct}</span>{middle_struct}"
        f"<span class='aug-span'>{aug_struct}</span>{after_struct}</div>"
    )


def make_downloads(dataset: Dict[str, Any]) -> None:
    candidates_df = to_df(dataset.get("allCandidateRows", []))
    binding_df = to_df(dataset.get("allBindingSiteRows", []))

    col1, col2, col3 = st.columns(3)

    with col1:
        st.download_button(
            "Download dataset JSON",
            data=json.dumps(dataset, indent=2),
            file_name="shine_dalguardyes_dataset.json",
            mime="application/json",
            use_container_width=True,
        )

    with col2:
        st.download_button(
            "Download candidates CSV",
            data=candidates_df.to_csv(index=False),
            file_name="shine_dalguardyes_candidates.csv",
            mime="text/csv",
            use_container_width=True,
            disabled=candidates_df.empty,
        )

    with col3:
        st.download_button(
            "Download binding sites CSV",
            data=binding_df.to_csv(index=False),
            file_name="shine_dalguardyes_binding_sites.csv",
            mime="text/csv",
            use_container_width=True,
            disabled=binding_df.empty,
        )


# -----------------------------------------------------------------------------
# Sidebar inputs
# -----------------------------------------------------------------------------

st.sidebar.markdown("# Shine DalguardYES")
st.sidebar.caption("ORBS-duplex to candidates.py to GA/TIR optimizer")

with st.sidebar.expander("Design inputs", expanded=True):
    orth_asd = st.text_input(
        "Orthogonal anti-SD",
        value="ACTTGTATA",
        help="This sequence is used as the ORBS input and as the orthogonal anti-SD.",
    )
    wt_asd = st.text_input("WT anti-SD", value="ACCTCCTTA")
    five_prime_flank = st.text_input("5′ flank / mRNA upstream", value="GCGGAAUUCGAUAA")
    cds_start = st.text_input("CDS start", value="AUGGCUACUAAAGAAAACGCU")
    job_name = st.text_input("Design name", value="shine_dalguardyes_design")

with st.sidebar.expander("ORBS-duplex settings", expanded=False):
    asd_tail_nt = st.number_input(
        "ASD tail nt to scan",
        min_value=4,
        max_value=50,
        value=12,
        step=1,
    )

    window_k = st.number_input(
        "Window k",
        min_value=3,
        max_value=12,
        value=4,
        step=1,
    )

    temp_c = st.number_input(
        "Temperature °C",
        min_value=20.0,
        max_value=45.0,
        value=37.0,
        step=0.5,
    )

with st.sidebar.expander("Candidate generator settings", expanded=False):
    candidates_per_seed = st.number_input(
        "Initial candidates from candidates.py",
        min_value=10,
        max_value=2000,
        value=120,
        step=10,
    )

    standby_start = st.number_input(
        "Frozen flank before index",
        min_value=0,
        max_value=200,
        value=0,
        step=1,
    )

    spacer_min, spacer_max = st.slider(
        "Spacer length range",
        min_value=3,
        max_value=14,
        value=(4, 7),
        step=1,
    )

    max_tries = st.number_input(
        "Max tries",
        min_value=1000,
        max_value=100000,
        value=10000,
        step=1000,
    )

with st.sidebar.expander("GA optimizer settings", expanded=True):
    population_size = st.number_input(
        "Population size",
        min_value=10,
        max_value=1000,
        value=80,
        step=10,
    )

    generations = st.number_input(
        "Generations",
        min_value=1,
        max_value=200,
        value=25,
        step=1,
    )

    elite_fraction = st.slider(
        "Elite fraction",
        min_value=0.05,
        max_value=0.50,
        value=0.20,
        step=0.05,
    )

    wt_penalty_constant = st.number_input(
        "WT penalty constant λ",
        min_value=0.0,
        max_value=10.0,
        value=1.5,
        step=0.1,
    )

    rng_seed = st.number_input(
        "Random seed",
        min_value=0,
        max_value=999999,
        value=7,
        step=1,
    )

    top_n = st.number_input(
        "Top candidates shown",
        min_value=5,
        max_value=100,
        value=20,
        step=5,
    )

run_clicked = st.sidebar.button(
    "Run Optimization",
    type="primary",
    use_container_width=True,
)


# -----------------------------------------------------------------------------
# Main dashboard
# -----------------------------------------------------------------------------

st.markdown(
    "<div class='hero-title'>Shine DalguardYES / Safety Dashboard</div>",
    unsafe_allow_html=True,
)

st.markdown(
    "<div class='hero-subtitle'>Orthogonal RBS Optimization and Biosafety Analysis</div>",
    unsafe_allow_html=True,
)
st.markdown(
    "<div class='hero-subtitle'>GitHub repo : <a href='https://github.com/crakshay1/DalguardYES' target='_blank'>crakshay1/DalguardYES</a></div>",
    unsafe_allow_html=True,
)

if not run_clicked and "dataset" not in st.session_state:
    st.info(
        "Set inputs in the sidebar, then click **Run Optimization**. "
        "The pipeline will run ORBS-duplex first, generate initial candidates, "
        "then optimize them with the GA/TIR engine."
    )

    vision_image = Path(__file__).parent / "synthetic_biology_optimization_dashboard_ui.png"

    if vision_image.exists():
        st.image(
            str(vision_image),
            caption="Target dashboard vision",
            use_container_width=True,
        )

    st.stop()

if run_clicked:
    try:
        with st.spinner("Running full pipeline: ORBS-duplex to candidates.py to GA/TIR optimizer..."):
            dataset = run_full_pipeline_cached(
                search_sequence=orth_asd,
                orth_asd=orth_asd,
                wt_asd=wt_asd,
                five_prime_flank=five_prime_flank,
                cds_start=cds_start,
                name=job_name,
                asd_tail_nt=int(asd_tail_nt),
                window_k=int(window_k),
                temp_c=float(temp_c),
                candidates_per_seed=int(candidates_per_seed),
                standby_start=int(standby_start),
                spacer_min=int(spacer_min),
                spacer_max=int(spacer_max),
                max_tries=int(max_tries),
                generations=int(generations),
                population_size=int(population_size),
                elite_fraction=float(elite_fraction),
                wt_penalty_constant=float(wt_penalty_constant),
                rng_seed=int(rng_seed),
                top_n=int(top_n),
            )

        st.session_state["dataset"] = dataset

    except Exception as exc:
        st.error(
            "Pipeline failed. Most common causes: ViennaRNA is not installed, "
            "ORBS search sequence is too short, or no candidates survived generation."
        )
        st.exception(exc)
        st.stop()


dataset = st.session_state["dataset"]
candidates = dataset.get("candidates", [])
all_candidates_df = to_df(dataset.get("allCandidateRows", []))
binding_sites_df = to_df(dataset.get("allBindingSiteRows", []))
fitness_df = to_df(dataset.get("fitnessData", []))
scatter_df = to_df(dataset.get("scatterPoints", []))

if not candidates or all_candidates_df.empty:
    st.warning("No candidates were returned.")
    st.stop()

best = candidates[0]
best_row = all_candidates_df.iloc[0].to_dict()
best_id = best_row["candidate_id"]


# -----------------------------------------------------------------------------
# KPI cards
# -----------------------------------------------------------------------------

k1, k2, k3, k4, k5, k6 = st.columns([1.15, 1, 1, 1, 1, 1])

with k1:
    st.markdown(
        f"""
        <div class='metric-card'>
          <div class='metric-label'>Best RBS</div>
          <div class='metric-value'>{best_row['rbs']}</div>
        </div>
        """,
        unsafe_allow_html=True,
    )

with k2:
    st.markdown(
        f"""
        <div class='metric-card'>
          <div class='metric-label'>Spacer Length</div>
          <div class='metric-value metric-value-blue'>{len(str(best_row['spacer']))} nt</div>
        </div>
        """,
        unsafe_allow_html=True,
    )

with k3:
    st.markdown(
        f"""
        <div class='metric-card'>
          <div class='metric-label'>Orthogonal TIR</div>
          <div class='metric-value'>{fmt_sci(best_row['orth_tir'])}</div>
        </div>
        """,
        unsafe_allow_html=True,
    )

with k4:
    st.markdown(
        f"""
        <div class='metric-card'>
          <div class='metric-label'>WT TIR / Leakage</div>
          <div class='metric-value metric-value-purple'>{fmt_sci(best_row['wt_tir'])}</div>
        </div>
        """,
        unsafe_allow_html=True,
    )

with k5:
    st.markdown(
        f"""
        <div class='metric-card'>
          <div class='metric-label'>T-score</div>
          <div class='metric-value metric-value-blue'>{fmt_sci(best_row['t_score'])}</div>
        </div>
        """,
        unsafe_allow_html=True,
    )

with k6:
    long_range = bool(best_row.get("long_range_flag", False))
    flag_text = "PASS" if not long_range else "FILTERED"

    st.markdown(
        f"""
        <div class='metric-card'>
          <div class='metric-label'>Structure Filter</div>
          <div class='metric-value'>{flag_text}</div>
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")


# -----------------------------------------------------------------------------
# Top dashboard row
# -----------------------------------------------------------------------------

left, middle, right = st.columns([1.05, 1.25, 0.9])

with left:
    st.markdown(
        "<div class='panel-title'>Pipeline Seed from ORBS-duplex</div>",
        unsafe_allow_html=True,
    )

    seed = dataset["pipeline"]["orbsSeed"]

    st.json(
        {
            #"full_rbs": seed_record.get("rbs"),
            "core": seed.get("core"),
            "flank": seed.get("five_prime_flank"),
            "cds_start": seed.get("cds_start"),
        },
        expanded=False,
    )

    with st.expander("View ORBS-duplex merged hits"):
        st.dataframe(
            to_df(dataset["pipeline"].get("orbsMergedHits", [])),
            use_container_width=True,
        )

with middle:
    st.markdown(
        "<div class='panel-title'>GA Evolution (T-score)</div>",
        unsafe_allow_html=True,
    )

    if not fitness_df.empty:
        fig = px.line(
            fitness_df,
            x="generation",
            y=["best", "avg"],
            labels={
                "value": "Normalized score",
                "generation": "Generation",
                "variable": "Trace",
            },
            height=310,
        )

        fig.update_layout(
            template="plotly_dark",
            paper_bgcolor="#0e1117",
            plot_bgcolor="#0e1117",
            margin=dict(l=10, r=10, t=10, b=10),
            legend=dict(orientation="h"),
        )

        st.plotly_chart(fig, use_container_width=True)

with right:
    st.markdown(
        "<div class='panel-title'>Energy Breakdown</div>",
        unsafe_allow_html=True,
    )

    energy_df = pd.DataFrame([
        {"term": "ΔG duplex", "value": best_row.get("dG_duplex_orth")},
        {"term": "ΔG start", "value": best_row.get("dG_start")},
        {"term": "ΔG standby", "value": best_row.get("dG_standby")},
        {"term": "ΔG spacing", "value": best_row.get("best_dG_spacing")},
        {"term": "ΔG mRNA unfolding", "value": best_row.get("best_dG_mrna_unfolding")},
        {"term": "ΔG total", "value": best_row.get("dG_total")},
    ])

    st.dataframe(
        energy_df,
        hide_index=True,
        use_container_width=True,
    )


# -----------------------------------------------------------------------------
# Candidate ranking and binding site inspector
# -----------------------------------------------------------------------------

st.markdown("---")

rank_col, site_col = st.columns([1.1, 1.25])

with rank_col:
    st.markdown(
        "<div class='panel-title'>Candidate Ranking</div>",
        unsafe_allow_html=True,
    )

    rank_cols = [
        "candidate_id",
        "rbs",
        "spacer",
        "orth_tir",
        "wt_tir",
        "t_score",
        "dG_total",
        "dG_duplex_orth",
        "best_dG_spacing",
        "best_dG_mrna_unfolding",
        "rbs_access",
        "long_range_flag",
    ]

    display_df = all_candidates_df[
        [c for c in rank_cols if c in all_candidates_df.columns]
    ].head(int(top_n)).copy()

    display_df.insert(0, "rank", range(1, len(display_df) + 1))

    st.dataframe(
        display_df,
        hide_index=True,
        use_container_width=True,
        height=360,
    )

    st.markdown(
        "<div class='small-note'>Expand rows below to inspect binding sites per candidate.</div>",
        unsafe_allow_html=True,
    )

    for _, row in all_candidates_df.head(min(int(top_n), 8)).iterrows():
        cid = row["candidate_id"]
        title = f"{cid} | RBS {row['rbs']} | T-score {fmt_sci(row['t_score'])}"

        with st.expander(title, expanded=(cid == best_id)):
            sub = binding_sites_df[
                (binding_sites_df["candidate_id"] == cid)
                & (binding_sites_df["anti_sd_type"] == "orthogonal")
            ].copy()

            keep_cols = [
                "site_rank",
                "mrna_start_1based",
                "asd_start_5p_1based",
                "aligned_spacing",
                "d",
                "dG_spacing",
                "dG_mrna_unfolding",
                "dG_total",
                "tir",
                "is_best_site",
            ]

            st.dataframe(
                sub[[c for c in keep_cols if c in sub.columns]],
                hide_index=True,
                use_container_width=True,
            )

with site_col:
    st.markdown(
        f"<div class='panel-title'>Best Candidate Binding Sites ({best_id})</div>",
        unsafe_allow_html=True,
    )

    best_sites = binding_sites_df[
        (binding_sites_df["candidate_id"] == best_id)
        & (binding_sites_df["anti_sd_type"] == "orthogonal")
    ].copy()

    site_cols = [
        "site_rank",
        "mrna_start_1based",
        "asd_start_5p_1based",
        "aligned_spacing",
        "d",
        "dG_spacing",
        "dG_mrna_unfolding",
        "dG_total",
        "tir",
        "is_best_site",
    ]

    st.dataframe(
        best_sites[[c for c in site_cols if c in best_sites.columns]],
        hide_index=True,
        use_container_width=True,
        height=255,
    )

    st.markdown(
        "<div class='panel-title'>Local Structure / Sequence Inspector</div>",
        unsafe_allow_html=True,
    )

    st.markdown(
        highlight_sequence(
            best_row["full_seq"],
            best_row["structure"],
            int(best_row["rbs_start"]),
            int(best_row["rbs_end"]),
            int(best_row["aug_start"]),
            int(best_row["aug_end"]),
        ),
        unsafe_allow_html=True,
    )


# -----------------------------------------------------------------------------
# Landscape and downloads
# -----------------------------------------------------------------------------

st.markdown("---")

land_col, dl_col = st.columns([1.25, 0.75])

with land_col:
    st.markdown(
        "<div class='panel-title'>Orthogonality Landscape</div>",
        unsafe_allow_html=True,
    )

    if not scatter_df.empty:
        fig = px.scatter(
            scatter_df,
            x="wtLeakage",
            y="binding",
            color="fitness",
            color_continuous_scale="Reds",
            hover_data=["rbs", "spacer"],
            labels={
                "wtLeakage": "WT TIR / leakage",
                "binding": "Orthogonal TIR",
                "fitness": "Fitness",
            },
            height=360,
            log_x=True,
            log_y=True,
        )

        fig.update_layout(
            template="plotly_dark",
            paper_bgcolor="#0e1117",
            plot_bgcolor="#0e1117",
            margin=dict(l=10, r=10, t=10, b=10),
        )

        st.plotly_chart(fig, use_container_width=True)

with dl_col:
    st.markdown(
        "<div class='panel-title'>Exports</div>",
        unsafe_allow_html=True,
    )

    make_downloads(dataset)

    st.markdown(
        "<div class='panel-title'>Pipeline Summary</div>",
        unsafe_allow_html=True,
    )

    st.write(" -> ".join(dataset["pipeline"]["steps"]))

    st.metric(
        "Initial candidates",
        len(dataset["pipeline"].get("initialCandidates", [])),
    )

    st.metric(
        "Ranked candidates",
        len(dataset.get("allCandidateRows", [])),
    )

    st.metric(
        "Binding-site rows",
        len(dataset.get("allBindingSiteRows", [])),
    )