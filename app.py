"""
RedoxScope - Polyphenol Toxicity, Docking & Apoptosis Explorer
Run with:  streamlit run app.py
"""
import datetime
import io
import zipfile

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

import helpers as H

APP_NAME = "RedoxScope"
APP_TAG = "Polyphenol Toxicity, Docking & Apoptosis Explorer"
PROTOX_URL = "https://tox.charite.de/protox3/"
CBDOCK_URL = "http://cao.labshare.cn:11080/cb-dock2/php/blinddock.php"

st.set_page_config(page_title=f"{APP_NAME} | {APP_TAG}", page_icon="🧪", layout="wide")

# ----------------------------------------------------------------------------
# Session state
# ----------------------------------------------------------------------------
ss = st.session_state
ss.setdefault("pairs", [])
ss.setdefault("form_n", 0)
ss.setdefault("submitted", False)
ss.setdefault("flash", None)
ss.setdefault("last_checks", None)

FIGS, TABLES = {}, {}  # collected during one run, used by the export tab
PCONFIG = {"displaylogo": False,
           "toImageButtonOptions": {"format": "png", "scale": 3, "filename": "redoxscope_plot"}}
ICON = {"ok": "✅", "warn": "🟡", "error": "❌"}


# ----------------------------------------------------------------------------
# UI helpers
# ----------------------------------------------------------------------------
def stretch(fn, *a, **k):
    """Works with old and new Streamlit versions (width='stretch' vs use_container_width)."""
    try:
        return fn(*a, width="stretch", **k)
    except Exception:  # noqa
        return fn(*a, use_container_width=True, **k)


def plot(fig, name):
    fig.update_layout(template="plotly_white", margin=dict(t=60, b=40, l=40, r=20))
    FIGS[name] = fig
    stretch(st.plotly_chart, fig, config=PCONFIG, key=f"plt_{name}")


def show_table(df, name, caption=None):
    TABLES[name] = df
    if caption:
        st.markdown(f"**{caption}**")
    stretch(st.dataframe, df, hide_index=True)
    st.download_button("⬇️ Download table (.csv)", df.to_csv(index=False).encode(),
                       file_name=f"{name}.csv", mime="text/csv", key=f"dl_{name}")


def summary_box(text):
    st.info("**📘 What is happening in this tab?**\n\n" + text)


def interpret_box(text):
    st.success("**🔎 How to read these results**\n\n" + text)


def show_checks(checks):
    for lvl, msg in checks:
        st.markdown(f"{ICON[lvl]} {msg}")


def gate():
    if not ss.submitted:
        st.warning("👈 Please finish **Tab 1** first: add at least **2 ligands** and press "
                   "**Submit & unlock analysis**.")
        return False
    return True


# ----------------------------------------------------------------------------
# Derived data
# ----------------------------------------------------------------------------
def build_data():
    pairs = ss.pairs
    D = {}
    ligs = {}
    for p in pairs:
        ligs.setdefault(p["ligand"], p)
    # ProTox endpoints (long + wide)
    long_frames = []
    for name, p in ligs.items():
        d = p["csv_df"].copy()
        d.insert(0, "Ligand", name)
        long_frames.append(d)
    D["tox_long"] = pd.concat(long_frames, ignore_index=True)
    sc = D["tox_long"].pivot_table(index="Ligand", columns="Shorthand",
                                   values="Active_Probability", aggfunc="first")
    for c in H.SHORTHANDS:
        if c not in sc.columns:
            sc[c] = np.nan
    D["scores"] = sc[H.SHORTHANDS].loc[list(ligs)]
    # PDF info
    D["info"] = pd.DataFrame([{
        "Ligand": n, "PDF compound name": p["pdf_info"].get("name"),
        "LD50 (mg/kg)": p["pdf_info"].get("ld50"), "Toxicity class": p["pdf_info"].get("tox_class"),
        "Mol. weight": p["pdf_info"].get("mol_weight"), "H-bond donors": p["pdf_info"].get("hbd"),
        "H-bond acceptors": p["pdf_info"].get("hba"), "logP": p["pdf_info"].get("logp"),
        "TPSA": p["pdf_info"].get("tpsa"), "Rotatable bonds": p["pdf_info"].get("rot_bonds")}
        for n, p in ligs.items()])
    # ProTox targets
    tf = []
    for n, p in ligs.items():
        t = p["targets_df"].copy()
        if not t.empty:
            t.insert(0, "Ligand", n)
            t["Pathway"] = t["Toxicity Target"].map(H.pathway_for)
            tf.append(t)
    D["targets"] = pd.concat(tf, ignore_index=True) if tf else pd.DataFrame()
    # Docking
    dk = []
    for p in pairs:
        d = p["dock_df"].copy()
        d.insert(0, "Protein", p["protein"])
        d.insert(0, "Ligand", p["ligand"])
        d.insert(0, "Pair", f"{p['ligand']} → {p['protein']}")
        dk.append(d)
    dock = pd.concat(dk, ignore_index=True)
    dock["Binding strength"] = dock["score"].map(H.strength_label)
    D["dock"] = dock
    best = dock.loc[dock.groupby("Pair")["score"].idxmin()].reset_index(drop=True)
    D["best"] = best
    return D


# ----------------------------------------------------------------------------
# File-stats bar (same on every tab)
# ----------------------------------------------------------------------------
def file_stats():
    pairs = ss.pairs
    ligs = {p["ligand"] for p in pairs}
    prots = {p["protein"] for p in pairs}
    cav = sorted({len(p["dock_df"]) for p in pairs})
    with st.container(border=True):
        c = st.columns(6)
        c[0].metric("Ligands", len(ligs))
        c[1].metric("Proteins", len(prots))
        c[2].metric("Ligand–protein pairs", len(pairs))
        c[3].metric("Files loaded", 3 * len(pairs))
        c[4].metric("Cavities per docking file", ", ".join(map(str, cav)) if cav else "–")
        c[5].metric("Status", "✅ Submitted" if ss.submitted else "⏳ Not submitted")
        if pairs:
            with st.expander("Show loaded files"):
                rows = [{"Ligand": p["ligand"], "Protein": p["protein"], "ProTox CSV": p["csv_name"],
                         "ProTox PDF": p["pdf_name"], "CurPockets TXT": p["dock_name"]} for p in pairs]
                stretch(st.dataframe, pd.DataFrame(rows), hide_index=True)


# ----------------------------------------------------------------------------
# Header + tabs
# ----------------------------------------------------------------------------
st.title(f"🧪 {APP_NAME}")
st.caption(f"**{APP_TAG}** · In-silico study: from antioxidant to pro-oxidant, and apoptosis prediction")

tabs = st.tabs(["1 · Upload & Check", "2 · Toxicity Comparison", "3 · Human Toxicity & Pathways",
                "4 · Docking Comparison", "5 · Antioxidant → Pro-oxidant", "6 · Apoptosis Prediction",
                "7 · Export Report"])


# ============================================================================
# TAB 1
# ============================================================================
def add_pair(ligand, pubchem, protein, pid, f_csv, f_pdf, f_txt):
    checks, ok_to_add = [], True
    missing = [n for n, v in [("ligand name", ligand), ("PubChem ID", pubchem), ("protein name", protein),
                              ("UniProt/PDB ID", pid), ("ProTox CSV", f_csv), ("ProTox PDF", f_pdf),
                              ("CurPockets TXT", f_txt)] if not v]
    if missing:
        return [("error", "Please provide: " + ", ".join(missing) + ".")], None

    csv_df, c1 = H.parse_protox_csv(f_csv.getvalue())
    info, targets, c2 = H.parse_protox_pdf(f_pdf.getvalue())
    dock, c3 = H.parse_curpockets(f_txt.getvalue())
    checks += [("ok" if l == "ok" else l, f"[CSV] {m}") for l, m in c1]
    checks += [(l, f"[PDF] {m}") for l, m in c2]
    checks += [(l, f"[TXT] {m}") for l, m in c3]
    if csv_df is None or info is None or dock is None:
        return checks, None

    # Ligand name vs name in the PDF
    pdf_name = info.get("name")
    if pdf_name:
        if H.name_similarity(ligand, pdf_name) < 0.5:
            checks.append(("warn", f"Ligand name '{ligand}' differs from the PDF compound name "
                                   f"'{pdf_name}'. Is this the right PDF?"))
        else:
            checks.append(("ok", f"Ligand name matches the ProTox PDF ('{pdf_name}')."))

    # Online ID mapping
    lvl, msg = H.check_pubchem(pubchem, ligand, info.get("mol_weight"))
    checks.append((lvl, msg))
    lvl, msg = H.check_protein_id(pid, protein)
    checks.append((lvl, msg))

    # Consistency with already added pairs
    for p in ss.pairs:
        if p["ligand"].lower() == ligand.lower() and p["protein"].lower() == protein.lower():
            checks.append(("error", "This ligand–protein pair is already in your list."))
        if p["ligand"].lower() == ligand.lower() and p["pubchem"] != pubchem.strip():
            checks.append(("error", f"'{ligand}' was already added with PubChem ID {p['pubchem']}."))
        if p["protein"].lower() == protein.lower() and p["protein_id"].upper() != pid.strip().upper():
            checks.append(("error", f"'{protein}' was already added with ID {p['protein_id']}."))
    if ss.pairs:
        expected = len(ss.pairs[0]["dock_df"])
        if len(dock) != expected:
            checks.append(("error", f"This docking file has {len(dock)} cavities but your earlier files have "
                                    f"{expected}. Use the SAME number of cavities for all files "
                                    f"(re-run CB-Dock2 if needed)."))
        else:
            checks.append(("ok", f"Cavity count ({len(dock)}) matches earlier files."))

    if any(l == "error" for l, _ in checks):
        return checks, None
    pair = {"ligand": ligand.strip(), "pubchem": pubchem.strip(), "protein": protein.strip(),
            "protein_id": pid.strip().upper(), "csv_df": csv_df, "csv_raw": f_csv.getvalue(),
            "csv_name": f_csv.name, "pdf_raw": f_pdf.getvalue(), "pdf_name": f_pdf.name,
            "pdf_info": info, "targets_df": targets, "dock_df": dock, "dock_raw": f_txt.getvalue(),
            "dock_name": f_txt.name,
            "mapped": all(l == "ok" for l, _ in checks)}
    return checks, pair


with tabs[0]:
    file_stats()
    st.header("Step 1 · Upload your files and check them")

    st.info(f"""
**📘 About {APP_NAME}.** This tool compares several *polyphenols* (plant compounds such as quercetin or genistein)
using two free web tools. For every **ligand (compound) + protein** pair you upload three files:

1. **ProTox-3 CSV** and 2. **ProTox-3 PDF** (toxicity prediction) → generate at **[{PROTOX_URL}]({PROTOX_URL})**
3. **CB-Dock2 `CurPockets_info.txt`** (docking result) → generate at **[{CBDOCK_URL}]({CBDOCK_URL})**
   and download the *CurPockets* file. {APP_NAME} works with the output of this website.

The app checks your files, converts the PDF toxicity-target table into a CSV, links your **PubChem ID** and
**UniProt/PDB ID** to the compounds and proteins, and then builds comparison plots in Tabs 2–6.
You can add as many ligand–protein pairs as you like, one by one (**at least 2 ligands are needed**).
""")
    with st.container(border=True):
        st.warning("""
**⚠️ Important notes before you start**
- **Use the same number of cavities for every CB-Dock2 analysis in one study.** A file with 5 cavities and
  another with 10 cavities will give unfair comparisons (more cavities = more chances of finding a good score).
  This app blocks files that do not match the first one you add.
- **Use the same ProTox-3 settings** for all ligands (same website, same date range/version, same options) so results are comparable.
- ProTox "Probability" is the *confidence of the predicted label*. This app converts it to **P(Active)** for plots.
- These are **computer predictions**, not lab results. They help you form hypotheses; they do not prove anything.
- File names can be anything, but each file must be the right *type* (CSV, PDF, TXT).
""")

    st.subheader("Add a ligand–protein pair")
    n = ss.form_n
    if ss.flash:
        st.success(ss.flash[0])
        show_checks(ss.flash[1])
        ss.flash = None
    c1, c2 = st.columns(2)
    with c1:
        st.markdown("##### 🧬 Ligand (compound)")
        ligand = st.text_input("Ligand name (e.g. Quercetin)", key=f"lig_{n}")
        pubchem = st.text_input("PubChem CID (e.g. 5280343)", key=f"cid_{n}")
        f_csv = st.file_uploader("ProTox-3 CSV file", type=["csv"], key=f"csv_{n}")
        f_pdf = st.file_uploader("ProTox-3 PDF file", type=["pdf"], key=f"pdf_{n}")
    with c2:
        st.markdown("##### 🧫 Protein used for docking")
        protein = st.text_input("Protein name (e.g. Prostaglandin G/H Synthase 1)", key=f"prot_{n}")
        pid = st.text_input("UniProt accession (e.g. P23219) or PDB ID", key=f"pid_{n}")
        f_txt = st.file_uploader("CB-Dock2 CurPockets_info TXT file", type=["txt"], key=f"txt_{n}")

    if st.button("🔍 Check files & add this pair", type="primary"):
        with st.spinner("Checking file formats and IDs..."):
            checks, pair = add_pair(ligand, pubchem, protein, pid, f_csv, f_pdf, f_txt)
        if pair:
            ss.pairs.append(pair)
            ss.submitted = False
            ss.form_n += 1
            ss.flash = (f"✅ Added **{pair['ligand']} → {pair['protein']}**. Add another pair or scroll down to submit.",
                        checks)
            st.rerun()
        else:
            st.error("Not added yet. Please fix the ❌ items below and try again.")
            show_checks(checks)

    # ---- Current list ----
    if ss.pairs:
        st.subheader("Your pairs so far")
        for i, p in enumerate(ss.pairs):
            cols = st.columns([5, 1])
            tick = "✅" if p["mapped"] else "🟡"
            cols[0].markdown(f"{tick} **{p['ligand']}** (CID {p['pubchem']}) → **{p['protein']}** "
                             f"({p['protein_id']}) · {len(p['dock_df'])} cavities")
            if cols[1].button("Remove", key=f"rm_{i}"):
                ss.pairs.pop(i)
                ss.submitted = False
                st.rerun()
        st.caption("✅ = formats and IDs confirmed · 🟡 = formats fine, but an online ID check could not be completed or looks different (please double-check).")

        # ---- Converted PDF tables ----
        st.subheader("PDF toxicity-target tables (converted to CSV)")
        for i, p in enumerate(dict((q["ligand"], q) for q in ss.pairs).values()):
            with st.expander(f"{p['ligand']}_Protox_toxicity targets.csv"):
                if p["targets_df"].empty:
                    st.write("No targets were listed in this PDF.")
                else:
                    stretch(st.dataframe, p["targets_df"], hide_index=True)
                    st.download_button("⬇️ Download CSV", p["targets_df"].to_csv(index=False).encode(),
                                       file_name=f"{H.safe_filename(p['ligand'])}_Protox_toxicity targets.csv",
                                       mime="text/csv", key=f"dlt_{i}")

        # ---- Submit ----
        st.subheader("Final check & submit")
        ligs = {p["ligand"] for p in ss.pairs}
        cav = {len(p["dock_df"]) for p in ss.pairs}
        problems = []
        if len(ligs) < 2:
            problems.append("Add at least **2 different ligands**.")
        if len(cav) > 1:
            problems.append("Docking files have different cavity counts.")
        for pr in problems:
            st.error(pr)
        if not ss.submitted:
            if st.button("🚀 Submit & unlock analysis", type="primary", disabled=bool(problems)):
                ss.submitted = True
                st.rerun()
        else:
            st.success("Submitted! Open Tabs 2–7 to explore the analysis.")
            if st.button("✏️ Edit data (unlock)"):
                ss.submitted = False
                st.rerun()

        D = build_data()
        best = D["best"].set_index("Pair")
        rows = []
        for i, p in enumerate(ss.pairs, 1):
            pair = f"{p['ligand']} → {p['protein']}"
            rows.append({"#": i, "Ligand": p["ligand"], "PubChem CID": p["pubchem"],
                         "Protein": p["protein"], "UniProt/PDB": p["protein_id"],
                         "ProTox CSV": f"{p['csv_name']} ({len(p['csv_df'])} rows)",
                         "ProTox PDF": f"{p['pdf_name']} ({len(p['targets_df'])} targets)",
                         "CurPockets TXT": f"{p['dock_name']} ({len(p['dock_df'])} cavities)",
                         "Best docking score": best.loc[pair, "score"],
                         "Mapping check": "✅" if p["mapped"] else "🟡"})
        show_table(pd.DataFrame(rows), "upload_summary", "Summary of all uploads")


# ============================================================================
# TAB 2: Toxicity comparison
# ============================================================================
def tab_toxicity(D):
    summary_box(
        "ProTox-3 predicts, with machine learning, whether each compound is **Active** (likely to cause an effect) "
        "or **Inactive** for several toxicity endpoints. Here we compare all ligands side by side. "
        "We use **P(Active)** (0 = surely inactive, 1 = surely active) so that all bars point the same way. "
        "We also show the predicted **LD50** (dose that kills 50% of rats: *higher = safer*) and a simple "
        "**drug-likeness** check (Lipinski's rule of five).")
    sc = D["scores"]
    wide = sc.copy()
    wide.columns = [H.NICE_NAMES[c] for c in wide.columns]
    wide = wide.round(2).reset_index()
    show_table(wide, "toxicity_P_active_table", "P(Active) for every endpoint")

    # Heatmap of toxicity endpoints
    hm = sc[H.TOX_ENDPOINTS].copy()
    hm.columns = [H.NICE_NAMES[c] for c in hm.columns]
    fig = px.imshow(hm, text_auto=".2f", color_continuous_scale="Reds", zmin=0, zmax=1, aspect="auto",
                    title="Toxicity endpoints heatmap: P(Active)", labels=dict(color="P(Active)"))
    plot(fig, "tox_heatmap")

    # Grouped bars
    long = hm.reset_index().melt(id_vars="Ligand", var_name="Endpoint", value_name="P(Active)")
    fig = px.bar(long, x="Endpoint", y="P(Active)", color="Ligand", barmode="group",
                 title="P(Active) by endpoint and ligand", range_y=[0, 1])
    plot(fig, "tox_grouped_bar")

    # Risk index
    risk = pd.DataFrame({"Ligand": sc.index, "Toxicity risk index": sc[H.TOX_ENDPOINTS].mean(axis=1).values,
                         "Endpoints flagged Active (P>0.5)": (sc[H.TOX_ENDPOINTS] > 0.5).sum(axis=1).values})
    risk = risk.sort_values("Toxicity risk index").reset_index(drop=True)
    risk["Toxicity risk index"] = risk["Toxicity risk index"].round(3)
    fig = px.bar(risk, x="Ligand", y="Toxicity risk index", text="Toxicity risk index", color="Toxicity risk index",
                 color_continuous_scale="OrRd", range_color=[0, 1], range_y=[0, 1],
                 title="Overall toxicity risk index (average P(Active) of 5 toxicity endpoints)")
    plot(fig, "tox_risk_index")
    show_table(risk, "toxicity_risk_index", "Ranking from lowest to highest risk")

    # LD50
    info = D["info"].copy()
    if info["LD50 (mg/kg)"].notna().any():
        info["Toxicity class"] = info["Toxicity class"].astype("Int64").astype(str)
        fig = px.bar(info, x="Ligand", y="LD50 (mg/kg)", color="Toxicity class", text="LD50 (mg/kg)",
                     title="Predicted oral LD50 (higher = less acutely toxic)",
                     color_discrete_sequence=px.colors.sequential.Viridis)
        plot(fig, "tox_ld50")

    # Lipinski
    li = D["info"].copy()
    li["MW ≤ 500"] = np.where(li["Mol. weight"] <= 500, "✅", "❌")
    li["logP ≤ 5"] = np.where(li["logP"] <= 5, "✅", "❌")
    li["Donors ≤ 5"] = np.where(li["H-bond donors"] <= 5, "✅", "❌")
    li["Acceptors ≤ 10"] = np.where(li["H-bond acceptors"] <= 10, "✅", "❌")
    li["Rules passed (of 4)"] = (li[["MW ≤ 500", "logP ≤ 5", "Donors ≤ 5", "Acceptors ≤ 10"]] == "✅").sum(axis=1)
    show_table(li, "druglikeness_lipinski", "Drug-likeness (Lipinski's rule of five) and compound properties")

    interpret_box(
        "- **Darker red = more likely to be toxic** for that endpoint. For a *cancer-treatment* idea you would like low "
        "liver (hepatotoxicity) and immune toxicity, so healthy cells are spared.\n"
        "- **Mutagenicity / carcinogenicity Active** means the compound might damage DNA. This can be a warning for safety, "
        "but DNA damage is also one way compounds push cancer cells toward death. Always discuss both sides.\n"
        "- **LD50:** a *lower* number means a *smaller* dose is dangerous. Toxicity class 1 is most toxic, class 6 is least toxic.\n"
        "- **Risk index** is a simple average made for this tool. Use it to *rank* ligands, not as a real safety number.\n"
        "- **Lipinski:** passing 4/4 rules suggests the compound could be absorbed when taken by mouth.")


# ============================================================================
# TAB 3: Human toxicity targets / pathways
# ============================================================================
def tab_pathways(D):
    summary_box(
        "ProTox-3 also checks whether your compound looks like known ligands of important human proteins "
        "(*toxicity targets*). We converted this table from your PDF. **Pharmacophore fit** tells how well the shape/chemical "
        "features of your compound fit the target; **similarity** tells how close your compound is to known binders. "
        "We link each target to a body pathway, then compare this with the **protein you chose for docking**.")
    T = D["targets"]
    if T.empty:
        st.info("No toxicity-target tables could be read from the PDFs.")
    else:
        T = T.copy()
        T["Concern level"] = pd.cut(T["Avg Pharmacophore Fit"], [-1, 20, 50, 101],
                                    labels=["Low (<20%)", "Medium (20–50%)", "High (>50%)"]).astype(str)
        show_table(T, "toxicity_targets_pathways", "Possible toxicity targets and affected pathways")

        sim = T.pivot_table(index="Toxicity Target", columns="Ligand", values="Avg Similarity Known Ligands")
        fig = px.imshow(sim, text_auto=".1f", color_continuous_scale="Oranges", aspect="auto", zmin=0, zmax=100,
                        title="Similarity to known ligands (%) per toxicity target", labels=dict(color="%"))
        plot(fig, "path_similarity_heatmap")

        fig = px.bar(T, x="Toxicity Target", y="Avg Pharmacophore Fit", color="Ligand", barmode="group",
                     title="Pharmacophore fit (%) per toxicity target")
        plot(fig, "path_fit_bar")

        cnt = T.groupby(["Pathway", "Ligand"]).size().reset_index(name="Number of targets")
        fig = px.bar(cnt, x="Pathway", y="Number of targets", color="Ligand", barmode="stack",
                     title="Which body pathways are possibly affected?")
        plot(fig, "path_pathway_counts")

        fig = px.scatter(T, x="Avg Similarity Known Ligands", y="Avg Pharmacophore Fit", color="Ligand",
                         symbol="Pathway", hover_name="Toxicity Target", size_max=14,
                         title="Fit vs similarity: top-right corner = strongest warning")
        fig.update_traces(marker=dict(size=13))
        plot(fig, "path_fit_vs_similarity")

    # Docked protein cross-check
    st.subheader("Docked protein vs ProTox toxicity targets")
    rows = []
    for _, b in D["best"].iterrows():
        row = {"Ligand": b["Ligand"], "Docked protein": b["Protein"], "Best cavity": int(b["ID"]),
               "Best docking score": b["score"], "Binding strength": b["Binding strength"],
               "In ProTox target list?": "No", "Matching ProTox target": "–", "Fit (%)": np.nan}
        if not T.empty:
            tl = T[T["Ligand"] == b["Ligand"]]
            if not tl.empty:
                m, s = H.best_match(b["Protein"], tl["Toxicity Target"].tolist())
                if s >= 0.6:
                    row["In ProTox target list?"] = "Yes ✅"
                    row["Matching ProTox target"] = m
                    row["Fit (%)"] = float(tl[tl["Toxicity Target"] == m]["Avg Pharmacophore Fit"].iloc[0])
        rows.append(row)
    cross = pd.DataFrame(rows)
    show_table(cross, "docked_protein_vs_protox", "Does the docked protein also appear as a ProTox toxicity target?")

    if D["best"]["Protein"].nunique() >= 1:
        hm = D["best"].pivot_table(index="Ligand", columns="Protein", values="score", aggfunc="min")
        fig = px.imshow(hm, text_auto=".1f", color_continuous_scale="Blues_r", aspect="auto",
                        title="Best docking score: ligands × proteins (more negative = stronger)")
        plot(fig, "path_ligand_protein_heatmap")

    interpret_box(
        "- **High fit (>50%) AND high similarity (>80%)** means the compound probably binds that protein and could disturb "
        "its pathway (side effects). Low fit values are less worrying.\n"
        "- **Same target for many ligands** (a column that is dark for every ligand) points to a *shared* effect of the polyphenol family.\n"
        "- If your **docked protein is also a ProTox target** and the docking score is strong (≤ −7), two independent methods "
        "agree: the compound likely really interacts with that protein.\n"
        "- If they disagree, that is not a mistake; the methods look at different things. Report both honestly.")


# ============================================================================
# TAB 4: Docking
# ============================================================================
def tab_docking(D):
    dock, best = D["dock"], D["best"]
    summary_box(
        "CB-Dock2 finds *cavities* (pockets) on each protein and docks your ligand into each. The **score** is an estimated "
        "binding energy in kcal/mol: **more negative = stronger binding**. Cavity **volume** is the pocket size (Å³). "
        "Here we compare the best pockets, all pockets, and which amino acids touch the ligand.")
    c = st.columns(3)
    top = best.loc[best["score"].idxmin()]
    c[0].metric("Strongest binder", f"{top['Ligand']}")
    c[1].metric("Best score", f"{top['score']} kcal/mol")
    c[2].metric("Docked protein", top["Protein"])

    show_table(dock.drop(columns=["contact_resi"]).round(2), "docking_all_cavities", "All cavities, all pairs")
    rank = best[["Pair", "ID", "Cavities_volume", "score", "Binding strength", "n_contacts"]].sort_values("score")
    rank.insert(0, "Rank", range(1, len(rank) + 1))
    rank = rank.rename(columns={"ID": "Best cavity", "n_contacts": "Contact residues"})
    show_table(rank, "docking_best_pocket_ranking", "Ranking of pairs by their best pocket")

    fig = px.bar(rank, x="Pair", y="score", color="Binding strength", text="score",
                 title="Best docking score per pair (lower = stronger)",
                 color_discrete_map={"Very strong": "#08519c", "Strong": "#3182bd", "Moderate": "#9ecae1", "Weak": "#d9d9d9"})
    fig.update_yaxes(autorange="reversed")
    plot(fig, "dock_best_bar")

    fig = px.bar(dock, x="ID", y="score", color="Pair", barmode="group",
                 title="Docking score in every cavity", labels={"ID": "Cavity ID", "score": "Score (kcal/mol)"})
    fig.update_xaxes(type="category")
    fig.update_yaxes(autorange="reversed")
    plot(fig, "dock_by_cavity")

    fig = px.scatter(dock, x="Cavities_volume", y="score", color="Pair", size="n_contacts", hover_data=["ID"],
                     title="Cavity volume vs docking score (bubble size = number of contact residues)",
                     labels={"Cavities_volume": "Cavity volume (Å³)", "score": "Score (kcal/mol)"})
    fig.update_yaxes(autorange="reversed")
    plot(fig, "dock_volume_vs_score")

    fig = px.box(dock, x="Pair", y="score", points="all", color="Pair", hover_data=["ID"],
                 title="Spread of scores across cavities")
    fig.update_yaxes(autorange="reversed")
    plot(fig, "dock_box")

    # Contact residues of best pocket
    st.subheader("Amino acids touching the ligand in the best pocket")
    sel = st.selectbox("Choose a pair", list(best["Pair"]), key="resi_pair")
    row = best[best["Pair"] == sel].iloc[0]
    res = pd.Series(H.residue_names(row["contact_resi"]))
    aa = res.str[:3].value_counts().reset_index()
    aa.columns = ["Amino acid", "Count"]
    fig = px.bar(aa, x="Amino acid", y="Count", title=f"Amino-acid types in the best pocket: {sel}")
    plot(fig, f"dock_residues_{H.safe_filename(sel)}")
    st.caption("Residues: " + ", ".join(res.tolist()))

    # Pocket overlap for same protein
    multi = best.groupby("Protein")["Ligand"].nunique()
    prots = multi[multi >= 2].index.tolist()
    if prots:
        st.subheader("Do different ligands use the same pocket?")
        rows = []
        for pr in prots:
            sub = best[best["Protein"] == pr]
            for i in range(len(sub)):
                for j in range(i + 1, len(sub)):
                    a, b = sub.iloc[i], sub.iloc[j]
                    ov = H.jaccard(H.residue_set(a["contact_resi"]), H.residue_set(b["contact_resi"]))
                    rows.append({"Protein": pr, "Ligand A": a["Ligand"], "Ligand B": b["Ligand"],
                                 "Best cavity A": int(a["ID"]), "Best cavity B": int(b["ID"]),
                                 "Shared residues (%)": round(100 * ov, 1)})
        ov = pd.DataFrame(rows)
        show_table(ov, "docking_pocket_overlap", "Overlap of contact residues in the best pockets")
        fig = px.bar(ov, x=ov["Ligand A"] + " vs " + ov["Ligand B"], y="Shared residues (%)", color="Protein",
                     title="Shared contact residues (%)", labels={"x": "Ligand pair"}, range_y=[0, 100])
        plot(fig, "dock_overlap_bar")

    interpret_box(
        "- Rough guide for scores: **≤ −9 very strong · −7 to −9 strong · −5 to −7 moderate · > −5 weak.** "
        "These are only rules of thumb.\n"
        "- The **best pocket** is the one with the lowest score; it is the most likely binding site.\n"
        "- If a ligand scores well in **many** cavities, it may stick to the protein in many places (less specific).\n"
        "- **Bigger pockets** can fit more contacts and often (not always) score better; check the scatter plot.\n"
        "- **High overlap (%)** between two ligands means they compete for the same pocket; low overlap means different binding sites.\n"
        "- Compare ligands **only on the same protein** and the **same number of cavities**.")


# ============================================================================
# TAB 5: Antioxidant -> pro-oxidant
# ============================================================================
def redox_table(D, thr):
    sc = D["scores"].fillna(0)
    out = pd.DataFrame({"Ligand": sc.index})
    out["Antioxidant defence (nrf2/ARE active)"] = sc["sr_are"].values
    out["Pro-oxidant stress index"] = sc[["sr_mmp", "sr_p53", "cyto"]].mean(axis=1).values
    out["Net redox balance (stress − defence)"] = out["Pro-oxidant stress index"] - out["Antioxidant defence (nrf2/ARE active)"]
    out["Redox profile"] = np.select(
        [out["Net redox balance (stress − defence)"] > thr, out["Net redox balance (stress − defence)"] < -thr],
        ["Pro-oxidant leaning", "Antioxidant leaning"], default="Balanced / transition zone")
    return out.round(3)


def tab_redox(D):
    summary_box(
        "Polyphenols are famous *antioxidants*, but at high concentrations some turn into **pro-oxidants**: they create "
        "reactive oxygen species (ROS) and stress the cell. ProTox-3 has no dose, so we use a simple teaching model:\n\n"
        "- **Antioxidant defence** = P(Active) of the **nrf2/ARE** pathway (the cell's own shield switching on).\n"
        "- **Pro-oxidant stress index** = average P(Active) of **MMP loss** (mitochondria stress), **p53** (DNA-damage alarm) and **cytotoxicity**.\n"
        "- **Net balance** = stress − defence. Positive = more stress than defence.\n\n"
        "*This is a transparent classroom model, not a validated clinical predictor.*")
    thr = st.slider("Transition-zone width (how big must the difference be to call it antioxidant or pro-oxidant?)",
                    0.05, 0.5, 0.2, 0.05, key="redox_thr")
    R = redox_table(D, thr)
    show_table(R, "redox_predictor_table", "Redox profile per ligand")

    long = R.melt(id_vars="Ligand", value_vars=["Antioxidant defence (nrf2/ARE active)", "Pro-oxidant stress index"],
                  var_name="Measure", value_name="Score")
    fig = px.bar(long, x="Ligand", y="Score", color="Measure", barmode="group", range_y=[0, 1],
                 title="Defence vs stress for each ligand")
    plot(fig, "redox_defence_vs_stress")

    fig = px.scatter(R, x="Antioxidant defence (nrf2/ARE active)", y="Pro-oxidant stress index", text="Ligand",
                     color="Redox profile", range_x=[-0.05, 1.05], range_y=[-0.05, 1.05],
                     title="Redox map: above the dotted line = more stress than defence")
    fig.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode="lines", line=dict(dash="dot", color="grey"),
                             showlegend=False, hoverinfo="skip"))
    fig.update_traces(textposition="top center", marker=dict(size=14), selector=dict(mode="markers+text"))
    plot(fig, "redox_map")

    fig = px.bar(R.sort_values("Net redox balance (stress − defence)"), x="Ligand",
                 y="Net redox balance (stress − defence)", color="Redox profile", text="Net redox balance (stress − defence)",
                 title="Net redox balance (positive = pro-oxidant leaning)",
                 color_discrete_map={"Pro-oxidant leaning": "#d7301f", "Antioxidant leaning": "#238b45",
                                     "Balanced / transition zone": "#fdae6b"})
    fig.add_hline(y=0, line_color="black")
    plot(fig, "redox_net_balance")

    interpret_box(
        "- **Bottom-right of the redox map** (high defence, low stress) → behaves like a classic antioxidant.\n"
        "- **Top-left** (low defence, high stress) → pro-oxidant leaning: the cell's shield is *not* switched on while "
        "mitochondria/p53/cytotoxicity signals are active.\n"
        "- **Near the diagonal** → balanced or in the *transition zone*, the exact place where dose decides whether a polyphenol helps or harms.\n"
        "- Example idea for a report: *'Ligand X shows low nrf2/ARE activity and high MMP activity, so at high doses it may act as a pro-oxidant.'*\n"
        "- Remember: ProTox is dose-independent. Real antioxidant-to-pro-oxidant switching depends on concentration, metal ions and cell type.")


# ============================================================================
# TAB 6: Apoptosis
# ============================================================================
def tab_apoptosis(D):
    summary_box(
        "**Apoptosis** is programmed cell death, the 'tidy' way cells self-destruct. A common path is: *stress/ROS ↑ → mitochondria lose "
        "membrane potential (MMP) → cytochrome c leaks out → caspases switch on → cell dies.* The tumour-suppressor **p53** can start "
        "this too. We combine four ProTox signals into a simple **Apoptosis Likelihood Score** (0–1). "
        "You can change the weights below.")
    with st.expander("⚙️ Change the weights (optional)"):
        c = st.columns(4)
        w_mmp = c[0].slider("MMP loss", 0, 100, 40, 5, key="w_mmp")
        w_p53 = c[1].slider("p53 activation", 0, 100, 20, 5, key="w_p53")
        w_cyto = c[2].slider("Cytotoxicity", 0, 100, 20, 5, key="w_cyto")
        w_are = c[3].slider("Shield failing (1 − nrf2/ARE)", 0, 100, 20, 5, key="w_are")
    w = np.array([w_mmp, w_p53, w_cyto, w_are], dtype=float)
    if w.sum() == 0:
        st.error("Please give at least one signal a weight above 0.")
        return
    w = w / w.sum()
    sc = D["scores"].fillna(0)
    sig = pd.DataFrame({"MMP loss": sc["sr_mmp"], "p53 activation": sc["sr_p53"], "Cytotoxicity": sc["cyto"],
                        "Shield failing (1 − nrf2/ARE)": 1 - sc["sr_are"]})
    contrib = sig * w
    A = pd.DataFrame({"Ligand": sig.index, "Apoptosis likelihood score": contrib.sum(axis=1).values})
    A["Category"] = pd.cut(A["Apoptosis likelihood score"], [-0.01, 0.4, 0.6, 1.01],
                           labels=["Low", "Moderate", "High"]).astype(str)
    A = pd.concat([A, sig.reset_index(drop=True)], axis=1).sort_values("Apoptosis likelihood score", ascending=False)
    A = A.round(3).reset_index(drop=True)
    A.insert(0, "Rank", range(1, len(A) + 1))
    show_table(A, "apoptosis_scores", "Apoptosis likelihood ranking")

    fig = px.bar(A, x="Ligand", y="Apoptosis likelihood score", color="Category", text="Apoptosis likelihood score",
                 range_y=[0, 1], title="Apoptosis Likelihood Score",
                 color_discrete_map={"Low": "#a1d99b", "Moderate": "#fdae6b", "High": "#d7301f"})
    plot(fig, "apop_score_bar")

    cl = (contrib.reset_index().melt(id_vars="Ligand", var_name="Signal", value_name="Contribution"))
    fig = px.bar(cl, x="Ligand", y="Contribution", color="Signal", barmode="stack",
                 title="What drives the score? (weighted contribution of each signal)")
    plot(fig, "apop_contributions")

    fig = go.Figure()
    cats = list(sig.columns)
    for lig in sig.index:
        vals = sig.loc[lig].tolist()
        fig.add_trace(go.Scatterpolar(r=vals + vals[:1], theta=cats + cats[:1], name=lig, fill="toself", opacity=0.5))
    fig.update_layout(polar=dict(radialaxis=dict(range=[0, 1])), title="Apoptosis signal fingerprint")
    plot(fig, "apop_radar")

    # Link to docking
    best_per_lig = D["best"].groupby("Ligand")["score"].min().reset_index(name="Best docking score")
    link = A.merge(best_per_lig, on="Ligand")
    link["Docking strength"] = -link["Best docking score"]
    fig = px.scatter(link, x="Docking strength", y="Apoptosis likelihood score", text="Ligand", color="Category",
                     range_y=[0, 1.05], title="Does the strongest binder also look the most apoptotic?",
                     labels={"Docking strength": "Docking strength (−score, higher = stronger binding)"},
                     color_discrete_map={"Low": "#a1d99b", "Moderate": "#fdae6b", "High": "#d7301f"})
    fig.update_traces(textposition="top center", marker=dict(size=14))
    plot(fig, "apop_vs_docking")
    show_table(link.drop(columns=["Docking strength"]), "apoptosis_vs_docking", "Apoptosis score with best docking score")

    interpret_box(
        "- **High (≥ 0.6)** means several signals point toward cell death: mitochondria stress, p53, cytotoxicity and a weak antioxidant shield.\n"
        "- **MMP loss (sr_mmp) is the most important signal** in this model: when mitochondria lose their electric charge the cell is pushed toward apoptosis.\n"
        "- The **stacked bar** shows *why* a ligand scored high; two ligands can reach the same score by different routes.\n"
        "- **Top-right of the last plot** = strong binding to your protein *and* high apoptosis score: a good candidate to discuss.\n"
        "- If cytotoxicity is high but p53 and MMP are low, cells may die by **necrosis** instead (messier, causes inflammation).\n"
        "- This is a ranking tool made for learning. Real apoptosis needs lab tests such as caspase or Annexin V assays.")
    return A


# ============================================================================
# TAB 7: Export
# ============================================================================
def build_report_html(findings, A, R):
    parts = [f"<html><head><meta charset='utf-8'><title>{APP_NAME} report</title>"
             "<script src='https://cdn.plot.ly/plotly-2.35.2.min.js'></script>"
             "<style>body{font-family:Arial;max-width:1000px;margin:auto;padding:20px}"
             "table{border-collapse:collapse}td,th{border:1px solid #ccc;padding:4px 8px}</style></head><body>",
             f"<h1>{APP_NAME} report</h1><p>{APP_TAG}<br>Generated {datetime.datetime.now():%Y-%m-%d %H:%M}</p>",
             "<h2>Key findings</h2><ul>" + "".join(f"<li>{f}</li>" for f in findings) + "</ul>"]
    for name, df in TABLES.items():
        parts.append(f"<h3>{name.replace('_', ' ')}</h3>" + df.to_html(index=False, border=0))
    for name, fig in FIGS.items():
        parts.append(fig.to_html(full_html=False, include_plotlyjs=False))
    parts.append("<p><i>Computer predictions only: ProTox-3 and CB-Dock2 results are not laboratory evidence.</i></p></body></html>")
    return "".join(parts)


def tab_export(D, A, R):
    best, dock = D["best"], D["dock"]
    top = best.loc[best["score"].idxmin()]
    n_lig, n_prot = D["scores"].shape[0], best["Protein"].nunique()
    risk = D["scores"][H.TOX_ENDPOINTS].mean(axis=1)
    safest, riskiest = risk.idxmin(), risk.idxmax()
    findings = [
        f"{n_lig} ligands and {n_prot} protein(s) were analysed ({len(ss.pairs)} ligand–protein pairs, "
        f"{len(dock)} docking results).",
        f"Strongest docking: <b>{top['Ligand']}</b> with <b>{top['Protein']}</b>, score {top['score']} kcal/mol "
        f"in cavity {int(top['ID'])} ({top['Binding strength'].lower()} binding).",
        f"Lowest overall toxicity-risk index: <b>{safest}</b> ({risk[safest]:.2f}); highest: <b>{riskiest}</b> ({risk[riskiest]:.2f}).",
    ]
    if A is not None and len(A):
        findings.append(f"Highest apoptosis likelihood score: <b>{A.iloc[0]['Ligand']}</b> ({A.iloc[0]['Apoptosis likelihood score']:.2f}, "
                        f"{A.iloc[0]['Category'].lower()}).")
    pro = R[R["Redox profile"] == "Pro-oxidant leaning"]["Ligand"].tolist()
    anti = R[R["Redox profile"] == "Antioxidant leaning"]["Ligand"].tolist()
    findings.append("Pro-oxidant leaning: " + (", ".join(pro) if pro else "none") +
                    " · Antioxidant leaning: " + (", ".join(anti) if anti else "none") + ".")
    if A is not None and len(A) and A.iloc[0]["Ligand"] == top["Ligand"]:
        findings.append(f"<b>{top['Ligand']}</b> is both the strongest binder and the highest apoptosis score, a consistent story for a report.")

    st.success("Analysis complete! Everything below is included in the download package.")
    c = st.columns(5)
    c[0].metric("Ligands", n_lig)
    c[1].metric("Proteins", n_prot)
    c[2].metric("Docking results", len(dock))
    c[3].metric("Best score", f"{top['score']}")
    c[4].metric("Plots in report", len(FIGS))

    st.subheader("📝 Overall results")
    for f in findings:
        st.markdown("- " + f.replace("<b>", "**").replace("</b>", "**"))

    summary_box(
        "The ZIP file contains every table as CSV, the converted *ligand_Protox_toxicity targets.csv* files, your original uploads, "
        "an interactive **HTML report** (all plots, open in any browser) and a short **README**. "
        "Plots can also be saved one by one as PNG using the camera icon on each plot's toolbar.")
    st.markdown("""
**How to write this up (simple template)**
1. *Question:* Do polyphenols switch from antioxidant to pro-oxidant and trigger apoptosis?
2. *Method:* ProTox-3 for toxicity, CB-Dock2 for docking (same cavity count for all), RedoxScope for scoring.
3. *Results:* Use the key findings above, plus 2–3 plots (Tabs 4, 5, 6).
4. *Limits:* computer predictions, no dose, no real cells. Say what you would test in the lab next.
""")

    # ZIP
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for i, (name, df) in enumerate(TABLES.items(), 1):
            z.writestr(f"tables/{i:02d}_{name}.csv", df.to_csv(index=False))
        done = set()
        for p in ss.pairs:
            if p["ligand"] not in done:
                done.add(p["ligand"])
                z.writestr(f"converted_pdf_tables/{H.safe_filename(p['ligand'])}_Protox_toxicity targets.csv",
                           p["targets_df"].to_csv(index=False))
            z.writestr(f"original_files/{H.safe_filename(p['ligand'])}__{H.safe_filename(p['protein'])}/{p['csv_name']}", p["csv_raw"])
            z.writestr(f"original_files/{H.safe_filename(p['ligand'])}__{H.safe_filename(p['protein'])}/{p['pdf_name']}", p["pdf_raw"])
            z.writestr(f"original_files/{H.safe_filename(p['ligand'])}__{H.safe_filename(p['protein'])}/{p['dock_name']}", p["dock_raw"])
        z.writestr("report/RedoxScope_report.html", build_report_html(findings, A, R))
        z.writestr("README.txt",
                   f"{APP_NAME} - {APP_TAG}\nGenerated: {datetime.datetime.now():%Y-%m-%d %H:%M}\n\n"
                   + "\n".join("- " + f.replace("<b>", "").replace("</b>", "") for f in findings)
                   + f"\n\nSources: ProTox-3 ({PROTOX_URL}) and CB-Dock2 ({CBDOCK_URL}).\n"
                     "All results are computer predictions, not laboratory evidence.\n")
    st.download_button("📦 Download complete analysis package (.zip)", buf.getvalue(),
                       file_name="RedoxScope_analysis.zip", mime="application/zip", type="primary")
    with st.expander("What is inside the ZIP?"):
        st.write(f"{len(TABLES)} tables · {len(ss.pairs)} sets of original files · converted PDF tables · interactive HTML report · README")


# ----------------------------------------------------------------------------
# Run tabs
# ----------------------------------------------------------------------------
D = build_data() if ss.pairs and ss.submitted else None
A = R = None
with tabs[1]:
    file_stats()
    st.header("Toxicity Comparison")
    if gate():
        tab_toxicity(D)
with tabs[2]:
    file_stats()
    st.header("Human Toxicity Targets & Pathways")
    if gate():
        tab_pathways(D)
with tabs[3]:
    file_stats()
    st.header("Docking Comparison")
    if gate():
        tab_docking(D)
with tabs[4]:
    file_stats()
    st.header("Antioxidant → Pro-oxidant Predictor")
    if gate():
        tab_redox(D)
        R = redox_table(D, ss.get("redox_thr", 0.2))
with tabs[5]:
    file_stats()
    st.header("Apoptosis Prediction")
    if gate():
        A = tab_apoptosis(D)
with tabs[6]:
    file_stats()
    st.header("Export Report")
    if gate():
        tab_export(D, A, R if R is not None else redox_table(D, 0.2))

st.divider()
st.caption(f"{APP_NAME} · for learning and hypothesis building only · ProTox-3 and CB-Dock2 are external tools by their respective authors.")
