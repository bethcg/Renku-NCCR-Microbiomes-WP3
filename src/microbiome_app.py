"""
Interactive dashboard for the WP3 wheat-rhizosphere microbiome analysis.

Renku serves this as a web-application session via the Procfile:
    web: streamlit run src/microbiome_app.py --server.port $RENKU_SESSION_PORT ...

Run locally with:   streamlit run src/microbiome_app.py
"""
import numpy as np
import plotly.graph_objects as go
import streamlit as st

import microbiome_data as md  # same folder (src/) is on sys.path when run from here

# ------------------------------------------------------------------ page setup
st.set_page_config(page_title="WP3 Microbiome Dashboard", layout="wide")

NAVY, GREEN, CYAN, MUTED = "#07182A", "#007057", "#02BACE", "#9FB3C8"

st.markdown(
    f"""
    <style>
      .stApp {{ background:{NAVY}; }}
      h1,h2,h3,h4,p,label,.stMarkdown {{ color:#FFFFFF !important; }}
      [data-testid="stMetricValue"] {{ color:{CYAN} !important; }}
      [data-testid="stMetricLabel"] {{ color:{MUTED} !important; }}
      section[data-testid="stSidebar"] {{ background:#0E2A40; }}
    </style>
    """,
    unsafe_allow_html=True,
)

st.title("Wheat rhizosphere microbiome — WP3")
st.markdown(
    "Sequential propagation of a *reproducible* wheat rhizosphere microbiome "
    "(Garrido-Sanz & Keel, UNIL). Data: "
    "[Zenodo 10.5281/zenodo.14514438](https://doi.org/10.5281/zenodo.14514438) "
    "(CC-BY-4.0). Served reproducibly by Renku."
)


@st.cache_data(show_spinner="Loading data…")
def get_data():
    return md.load()


counts_all, meta, tax = get_data()
order = md.group_order(meta)
pal = md.palette(meta)

# ------------------------------------------------------------------ sidebar
st.sidebar.header("Controls")
groups = st.sidebar.multiselect("Community stages", order, default=order)
metric = st.sidebar.radio("Alpha-diversity metric", ["Shannon", "Observed"], horizontal=True)
min_prev = st.sidebar.slider("Min. prevalence (samples)", 1, 10, 2)
min_total = st.sidebar.slider("Min. total count per ASV", 0, 100, 10, step=5)
show_traj = st.sidebar.toggle("Show assembly trajectory", value=True)

if not groups:
    st.warning("Select at least one community stage in the sidebar.")
    st.stop()

# subset + filter
keep_samples = meta.index[meta["Name"].isin(groups)]
counts = md.filter_counts(counts_all.loc[keep_samples], min_prev, min_total)
meta_s = meta.loc[keep_samples]
ordered_groups = [g for g in order if g in groups]

# ------------------------------------------------------------------ KPIs
c1, c2, c3, c4 = st.columns(4)
c1.metric("Samples", f"{counts.shape[0]}")
c2.metric("ASVs (after filter)", f"{counts.shape[1]:,}")
c3.metric("Community stages", f"{len(ordered_groups)}")
c4.metric("Replicates / stage", f"{int(meta_s.groupby('Name').size().median())}")

alpha = md.alpha_table(counts, meta_s)
ord_df, var = md.braycurtis_pcoa(counts, meta_s)

# ------------------------------------------------------------------ figures
left, right = st.columns(2)

with left:
    st.subheader("Alpha diversity")
    fig = go.Figure()
    for g in ordered_groups:
        sub = alpha[alpha["Name"] == g]
        fig.add_trace(go.Box(y=sub[metric], name=g, marker_color=pal[g],
                             boxpoints="all", jitter=0.4, pointpos=0,
                             line_color=pal[g], fillcolor=pal[g], opacity=0.85))
    fig.update_layout(
        template="plotly_dark", paper_bgcolor=NAVY, plot_bgcolor=NAVY,
        showlegend=False, height=460, margin=dict(l=10, r=10, t=30, b=10),
        yaxis_title=f"{metric} diversity",
        xaxis=dict(categoryorder="array", categoryarray=ordered_groups))
    st.plotly_chart(fig, use_container_width=True)

with right:
    st.subheader(f"Bray-Curtis PCoA · PCo1 {var[0]:.1f}% · PCo2 {var[1]:.1f}%")
    fig = go.Figure()
    if show_traj and len(ordered_groups) > 1:
        cen = ord_df.groupby("Name")[["PCo1", "PCo2"]].mean().reindex(ordered_groups)
        fig.add_trace(go.Scatter(x=cen["PCo1"], y=cen["PCo2"], mode="lines",
                                 line=dict(color=MUTED, dash="dash"),
                                 hoverinfo="skip", showlegend=False))
    for g in ordered_groups:
        sub = ord_df[ord_df["Name"] == g]
        fig.add_trace(go.Scatter(
            x=sub["PCo1"], y=sub["PCo2"], mode="markers", name=g,
            marker=dict(size=11, color=pal[g], line=dict(width=0.5, color="black")),
            text=sub["id_samples"],
            hovertemplate="<b>%{text}</b><br>" + g + "<extra></extra>"))
    fig.update_layout(
        template="plotly_dark", paper_bgcolor=NAVY, plot_bgcolor=NAVY,
        height=460, margin=dict(l=10, r=10, t=30, b=10),
        legend=dict(font=dict(size=10)),
        xaxis_title=f"PCo1 ({var[0]:.1f}%)", yaxis_title=f"PCo2 ({var[1]:.1f}%)")
    st.plotly_chart(fig, use_container_width=True)

st.caption(
    "Left: diversity collapses through the seed-borne bottleneck. "
    "Right: the four replicates of each stage cluster tightly (reproducible) "
    "and trace the community-assembly trajectory. Same result every time it is "
    "re-run — that is the reproducibility Renku guarantees."
)

with st.expander("Per-sample diversity table"):
    st.dataframe(alpha[["id_samples", "Name", "Observed", "Shannon"]]
                 .reset_index(drop=True), use_container_width=True)
