"""RedoxScope helpers: file parsing, format checks, ID checks and simple scoring.
No Streamlit code in here, so everything can be tested on its own."""
import difflib
import io
import re

import numpy as np
import pandas as pd
import requests

# ----------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------
CSV_COLS = ["Classification", "Target", "Shorthand", "Prediction", "Probability"]
SHORTHANDS = ["dili", "carcino", "immuno", "mutagen", "cyto", "sr_are", "sr_mmp", "sr_p53"]
TOX_ENDPOINTS = ["dili", "carcino", "immuno", "mutagen", "cyto"]
NICE_NAMES = {
    "dili": "Hepatotoxicity (liver)",
    "carcino": "Carcinogenicity",
    "immuno": "Immunotoxicity",
    "mutagen": "Mutagenicity",
    "cyto": "Cytotoxicity",
    "sr_are": "nrf2/ARE antioxidant response",
    "sr_mmp": "Mitochondrial membrane potential (MMP)",
    "sr_p53": "p53 tumour suppressor",
}
POCKET_COLS = ["ID", "Cavities_volume", "center_x", "center_y", "center_z",
               "size_x", "size_y", "size_z", "score", "contact_resi"]

# Simple, beginner-friendly mapping of ProTox toxicity targets to biological pathways
PATHWAY_MAP = {
    "androgen receptor": "Hormone signalling",
    "estrogen receptor": "Hormone signalling",
    "progesterone receptor": "Hormone signalling",
    "glucocorticoid receptor": "Hormone signalling",
    "pregnane x receptor": "Drug metabolism / detox",
    "amine oxidase": "Neurotransmitter breakdown (brain)",
    "prostaglandin": "Inflammation & pain",
    "adenosine": "Nerve & heart signalling",
    "adrenergic": "Nerve & heart signalling",
    "dopamine": "Nerve & heart signalling",
    "histamine": "Allergy & immune response",
    "opioid": "Pain & nervous system",
    "corticotropin": "Stress response",
    "phosphodiesterase": "Cell signalling (cAMP)",
}


# ----------------------------------------------------------------------------
# Small utilities
# ----------------------------------------------------------------------------
def pathway_for(target_name: str) -> str:
    t = target_name.lower()
    for key, val in PATHWAY_MAP.items():
        if key in t:
            return val
    return "Other / not classified"


def _tokens(s):
    return {t for t in re.findall(r"[a-z0-9]+", str(s).lower()) if len(t) > 2}


def name_similarity(a, b):
    ta, tb = _tokens(a), _tokens(b)
    if not ta or not tb:
        return 0.0
    jac = len(ta & tb) / len(ta | tb)
    seq = difflib.SequenceMatcher(None, str(a).lower(), str(b).lower()).ratio()
    return max(jac, seq)


def best_match(name, candidates):
    """Return (best_candidate, similarity) for a name among candidates."""
    best, score = None, 0.0
    for c in candidates:
        s = name_similarity(name, c)
        if s > score:
            best, score = c, s
    return best, score


def safe_filename(s):
    return re.sub(r"[^A-Za-z0-9_.\- ]+", "_", str(s)).strip() or "ligand"


def strength_label(score):
    """Rule-of-thumb labels for AutoDock-Vina-style scores (kcal/mol)."""
    if score <= -9:
        return "Very strong"
    if score <= -7:
        return "Strong"
    if score <= -5:
        return "Moderate"
    return "Weak"


def tox_class_from_ld50(ld50):
    if ld50 <= 5:
        return 1
    if ld50 <= 50:
        return 2
    if ld50 <= 300:
        return 3
    if ld50 <= 2000:
        return 4
    if ld50 <= 5000:
        return 5
    return 6


# ----------------------------------------------------------------------------
# File parsers. Each returns (data, issues) where issues = [(level, message)]
# level is "ok", "warn" or "error"
# ----------------------------------------------------------------------------
def parse_protox_csv(raw: bytes):
    issues = []
    try:
        df = pd.read_csv(io.BytesIO(raw))
    except Exception as e:  # noqa
        return None, [("error", f"ProTox CSV could not be read: {e}")]
    df.columns = [str(c).strip() for c in df.columns]
    missing = [c for c in CSV_COLS if c not in df.columns]
    if missing:
        return None, [("error", f"ProTox CSV is missing column(s): {', '.join(missing)}. "
                                f"Expected: {', '.join(CSV_COLS)}")]
    df = df[CSV_COLS].copy()
    for c in ["Classification", "Target", "Shorthand", "Prediction"]:
        df[c] = df[c].astype(str).str.strip()
    df["Probability"] = pd.to_numeric(df["Probability"], errors="coerce")
    if df["Probability"].isna().any():
        return None, [("error", "ProTox CSV has non-numeric values in the Probability column.")]
    bad_pred = ~df["Prediction"].isin(["Active", "Inactive"])
    if bad_pred.any():
        return None, [("error", "ProTox CSV 'Prediction' must be Active or Inactive.")]
    missing_sh = [s for s in SHORTHANDS if s not in set(df["Shorthand"])]
    if missing_sh:
        issues.append(("warn", f"Some endpoints are missing in the CSV ({', '.join(missing_sh)}). "
                               "Related plots will treat them as 0."))
    # ProTox 'Probability' = confidence in the predicted label. Convert to P(Active).
    df["Active_Probability"] = np.where(df["Prediction"] == "Active",
                                        df["Probability"], 1 - df["Probability"])
    issues.insert(0, ("ok", f"ProTox CSV format OK ({len(df)} endpoints)."))
    return df, issues


def _first(pattern, text, flags=0, cast=float):
    m = re.search(pattern, text, flags)
    if not m:
        return None
    try:
        return cast(m.group(1))
    except Exception:  # noqa
        return m.group(1).strip()


def parse_protox_pdf(raw: bytes):
    """Returns (info_dict, targets_df, issues)."""
    import pdfplumber

    issues = []
    try:
        with pdfplumber.open(io.BytesIO(raw)) as pdf:
            text = "\n".join((p.extract_text() or "") for p in pdf.pages)
    except Exception as e:  # noqa
        return None, None, [("error", f"ProTox PDF could not be read: {e}")]

    if "ProTox" not in text and "Predicted LD50" not in text:
        return None, None, [("error", "This PDF does not look like a ProTox-3 result page.")]

    info = {
        "name": _first(r"\bName\s+([A-Za-z0-9][^\n]*)", text, cast=str),
        "ld50": _first(r"Predicted LD50:\s*([\d.]+)", text),
        "tox_class": _first(r"Predicted Toxicity Class:\s*(\d)", text, cast=int),
        "mol_weight": _first(r"Molweight\s+([\d.]+)", text),
        "hba": _first(r"Number of hydrogen\s+(\d+)\s*\n\s*bond acceptors", text, cast=int),
        "hbd": _first(r"Number of hydrogen\s+(\d+)\s*\n\s*bond donors", text, cast=int),
        "rot_bonds": _first(r"rotable bonds\s+(\d+)", text, cast=int),
        "tpsa": _first(r"Topological Polar\s+([\d.]+)", text),
        "logp": _first(r"octanol/water partition\s+(-?[\d.]+)", text),
    }
    if info["ld50"] is None or info["tox_class"] is None:
        issues.append(("warn", "Could not read LD50 / toxicity class from the PDF."))

    # --- Table: 'Details about possible toxicity targets' ---
    rows, capture = [], False
    row_re = re.compile(r"^(.*?)\s+(\d+(?:\.\d+)?)%\s+(\d+(?:\.\d+)?)%\s*$")
    for line in text.split("\n"):
        line = line.strip()
        if "Toxicity Target" in line and "Avg Pharmacophore Fit" in line:
            capture = True
            continue
        if capture:
            m = row_re.match(line)
            if m:
                rows.append({
                    "Toxicity Target": m.group(1).strip(),
                    "Avg Pharmacophore Fit": float(m.group(2)),
                    "Avg Similarity Known Ligands": float(m.group(3)),
                })
    targets = pd.DataFrame(rows, columns=["Toxicity Target", "Avg Pharmacophore Fit",
                                          "Avg Similarity Known Ligands"])
    if targets.empty:
        issues.append(("warn", "No rows found in 'Details about possible toxicity targets'. "
                               "(ProTox sometimes lists none.) Tab 3 will skip this ligand."))
    else:
        issues.insert(0, ("ok", f"ProTox PDF read OK; toxicity-target table converted "
                                f"({len(targets)} targets)."))
    return info, targets, issues


def parse_curpockets(raw: bytes):
    issues = []
    try:
        df = pd.read_csv(io.BytesIO(raw), sep=r"\s+", engine="python")
    except Exception as e:  # noqa
        return None, [("error", f"CurPockets file could not be read: {e}")]
    df.columns = [str(c).strip() for c in df.columns]
    missing = [c for c in POCKET_COLS if c not in df.columns]
    if missing:
        return None, [("error", f"CurPockets file is missing column(s): {', '.join(missing)}. "
                                "Please upload the *CurPockets_info.txt* file from CB-Dock2.")]
    if df.empty:
        return None, [("error", "CurPockets file has no cavities (rows).")]
    for c in POCKET_COLS[:-1]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    if df[POCKET_COLS[:-1]].isna().any().any():
        return None, [("error", "CurPockets file has non-numeric values in numeric columns.")]
    df["contact_resi"] = df["contact_resi"].astype(str)
    df["n_contacts"] = df["contact_resi"].str.count(",") + 1
    issues.append(("ok", f"CurPockets format OK ({len(df)} cavities)."))
    return df, issues


def residue_set(contact_str):
    return set(str(contact_str).split(","))


def residue_names(contact_str):
    """'TYR:38:A,ALA:89:A' -> ['TYR38','ALA89']"""
    out = []
    for r in str(contact_str).split(","):
        p = r.split(":")
        if len(p) >= 2:
            out.append(f"{p[0]}{p[1]}")
    return out


# ----------------------------------------------------------------------------
# Online ID checks (PubChem, UniProt, RCSB). They never crash the app.
# ----------------------------------------------------------------------------
UNIPROT_RE = re.compile(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$")
PDB_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")


def check_pubchem(cid: str, ligand_name: str, pdf_mw):
    cid = str(cid).strip()
    if not cid.isdigit():
        return "error", "PubChem ID must be a number (e.g. 5280343)."
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/Title,MolecularWeight/JSON"
    try:
        r = requests.get(url, timeout=10)
    except Exception:  # noqa
        return "warn", "PubChem ID format OK, but PubChem could not be reached to confirm it."
    if r.status_code == 404:
        return "error", f"PubChem CID {cid} was not found."
    if r.status_code != 200:
        return "warn", "PubChem ID format OK, but PubChem did not answer; mapping not confirmed."
    try:
        p = r.json()["PropertyTable"]["Properties"][0]
        title, mw = p.get("Title", ""), float(p.get("MolecularWeight"))
    except Exception:  # noqa
        return "warn", "PubChem answered in an unexpected way; mapping not confirmed."
    if pdf_mw is not None and abs(mw - pdf_mw) > 1.0:
        return "error", (f"PubChem CID {cid} is '{title}' (MW {mw:.2f}) but the ProTox PDF compound has "
                         f"MW {pdf_mw:.2f}. Wrong PubChem ID or wrong file?")
    sim = name_similarity(ligand_name, title)
    if sim < 0.5:
        return "warn", (f"PubChem CID {cid} = '{title}', MW matches the PDF, but the name differs from "
                        f"'{ligand_name}'. Please double-check (synonyms are fine).")
    return "ok", f"PubChem CID {cid} = {title} (MW matches the ProTox PDF)."


def check_protein_id(pid: str, protein_name: str):
    pid = str(pid).strip().upper()
    if UNIPROT_RE.match(pid):
        url = f"https://rest.uniprot.org/uniprotkb/{pid}?fields=protein_name,gene_names,organism_name&format=json"
        try:
            r = requests.get(url, timeout=10)
        except Exception:  # noqa
            return "warn", "UniProt ID format OK, but UniProt could not be reached to confirm it."
        if r.status_code in (400, 404):
            return "error", f"UniProt accession {pid} was not found."
        if r.status_code != 200:
            return "warn", "UniProt ID format OK, but UniProt did not answer; mapping not confirmed."
        try:
            j = r.json()
            pdsc = j.get("proteinDescription", {})
            names = []
            for n in ([pdsc.get("recommendedName")] if pdsc.get("recommendedName") else []) + \
                    pdsc.get("alternativeNames", []) + pdsc.get("submissionNames", []):
                v = n.get("fullName", {}).get("value")
                if v:
                    names.append(v)
            for g in j.get("genes", []):
                v = g.get("geneName", {}).get("value")
                if v:
                    names.append(v)
            org = j.get("organism", {}).get("scientificName", "")
        except Exception:  # noqa
            return "warn", "UniProt answered in an unexpected way; mapping not confirmed."
        _, sim = best_match(protein_name, names)
        main = names[0] if names else "unknown"
        if sim < 0.5:
            return "warn", (f"UniProt {pid} = '{main}' ({org}). This does not obviously match "
                            f"'{protein_name}'. Please double-check.")
        return "ok", f"UniProt {pid} = {main} ({org}); matches the protein name."
    if PDB_RE.match(pid):
        url = f"https://data.rcsb.org/rest/v1/core/entry/{pid.lower()}"
        try:
            r = requests.get(url, timeout=10)
        except Exception:  # noqa
            return "warn", "PDB ID format OK, but RCSB could not be reached to confirm it."
        if r.status_code == 404:
            return "error", f"PDB ID {pid} was not found at RCSB."
        if r.status_code != 200:
            return "warn", "PDB ID format OK, but RCSB did not answer; mapping not confirmed."
        try:
            title = r.json()["struct"]["title"]
        except Exception:  # noqa
            return "warn", "RCSB answered in an unexpected way; mapping not confirmed."
        if name_similarity(protein_name, title) < 0.3:
            return "warn", f"PDB {pid} = '{title}'. Please check this is your protein."
        return "ok", f"PDB {pid} = {title}."
    return "error", "ID is neither a valid UniProt accession (e.g. P23219) nor a 4-character PDB ID (e.g. 6Y3C)."


# ----------------------------------------------------------------------------
# Simple teaching scores (explained to the user in the app)
# ----------------------------------------------------------------------------
def jaccard(a, b):
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if (a | b) else 0.0
