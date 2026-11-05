# 🧪 RedoxScope

**Polyphenol Toxicity, Docking & Apoptosis Explorer**

RedoxScope is a beginner-friendly [Streamlit](https://streamlit.io) app that compares several polyphenols (plant compounds such as quercetin or genistein) using the output of two free web tools:

| Tool | Used for | Link |
|---|---|---|
| **ProTox-3** | Toxicity prediction (CSV + PDF) | https://tox.charite.de/protox3/ |
| **CB-Dock2** | Blind docking (`CurPockets_info.txt`) | http://cao.labshare.cn:11080/cb-dock2/php/blinddock.php |

It checks your files, links your PubChem and UniProt/PDB IDs to the compounds and proteins, and builds interactive comparison plots about toxicity, protein binding, antioxidant vs pro-oxidant behaviour, and apoptosis (programmed cell death) likelihood. It is designed so that high school students can use it.

> ⚠️ All results are **computer predictions** for learning and hypothesis building. They are not laboratory evidence.

---

## What you need

For **each ligand–protein pair** you need three files:

1. **ProTox-3 CSV** (the "CSV" button on the ProTox-3 results page)
2. **ProTox-3 PDF** (print/save the results page as PDF)
3. **CB-Dock2 CurPockets file**, named like `Quercetin_Protein_CurPockets_info.txt`

You also need:
- Ligand name and **PubChem CID** (e.g. Quercetin = `5280343`)
- Protein name and **UniProt accession** (e.g. `P23219`) or **PDB ID**

**At least 2 ligands** are required. A ligand can be docked to several proteins, and a protein can be used for several ligands. File names can be anything.

---

## ⚠️ Rules for a fair comparison

- **Use the same number of cavities in every CB-Dock2 analysis.** A file with 5 cavities and one with 10 cannot be compared fairly. The app rejects files that do not match the first one you add.
- **Use the same ProTox-3 settings** for all ligands.
- ProTox "Probability" is the confidence in the predicted label. RedoxScope converts it to **P(Active)** before plotting.

---

## The 7 tabs

| Tab | What it does |
|---|---|
| **1 · Upload & Check** | Add ligand–protein pairs one by one. Checks file formats, converts the PDF "Details about possible toxicity targets" table to `<ligand>_Protox_toxicity targets.csv`, confirms PubChem and UniProt/PDB IDs online (✅ / 🟡 / ❌), then a final summary table and **Submit**. |
| **2 · Toxicity Comparison** | Heatmap and bars of toxicity endpoints, overall risk index, predicted LD50 and toxicity class, drug-likeness (Lipinski) table. |
| **3 · Human Toxicity & Pathways** | Toxicity targets from the PDF mapped to body pathways; checks whether your docked protein is also a ProTox target. |
| **4 · Docking Comparison** | Best-pocket ranking, scores per cavity, volume vs score, score spread, contact amino acids, and pocket overlap between ligands on the same protein. |
| **5 · Antioxidant → Pro-oxidant** | Defence (nrf2/ARE) vs stress (MMP, p53, cytotoxicity) map and net redox balance with an adjustable threshold. |
| **6 · Apoptosis Prediction** | Apoptosis Likelihood Score (adjustable weights), contribution chart, radar plot, link to docking strength. |
| **7 · Export Report** | Overall summary, key findings, write-up template, and a **ZIP** with all tables, converted PDF tables, original files, an interactive HTML report and a README. |

Every tab shows the same file-stats bar, a "What is happening" block, a "How to read these results" block, Plotly charts (with toolbar and PNG download) and CSV downloads for each table.

---

## Run locally

```bash
git clone https://github.com/<your-username>/redoxscope.git
cd redoxscope
pip install -r requirements.txt
streamlit run app.py
```

Python 3.9+ is recommended. Online ID checks need an internet connection; without it the app still works and marks IDs as 🟡 (format OK, not confirmed).

## Deploy on Streamlit Community Cloud

1. Push `app.py`, `helpers.py`, `requirements.txt` and `README.md` to a GitHub repository.
2. Go to [share.streamlit.io](https://share.streamlit.io), choose **New app**, select the repository, and set the main file to `app.py`.
3. Click **Deploy**.

## Project files

```
app.py            Streamlit interface and all 7 tabs
helpers.py        File parsing, format checks, ID checks, simple scoring (no Streamlit code)
requirements.txt  Python dependencies
```

---

## How the simple scores work

These are transparent teaching models, **not validated predictors**.

- **P(Active)** = Probability if Prediction is *Active*, otherwise 1 − Probability.
- **Toxicity risk index** = mean P(Active) of hepatotoxicity, carcinogenicity, immunotoxicity, mutagenicity and cytotoxicity.
- **Antioxidant defence** = P(Active) of nrf2/ARE.
- **Pro-oxidant stress index** = mean P(Active) of MMP, p53 and cytotoxicity.
- **Net redox balance** = stress − defence (positive = pro-oxidant leaning; the transition zone width is set by a slider).
- **Apoptosis Likelihood Score** = weighted mean of MMP (40%), p53 (20%), cytotoxicity (20%) and "shield failing" = 1 − nrf2/ARE (20%). Weights can be changed in Tab 6.
- **Docking strength guide** (kcal/mol): ≤ −9 very strong, −7 to −9 strong, −5 to −7 moderate, > −5 weak.

ProTox-3 gives no dose information, so the antioxidant-to-pro-oxidant switch is only a hint. Real switching depends on concentration, metal ions and cell type.

---

## Limitations

- Predictions only; no dose, no real cells. Confirm ideas with laboratory work (e.g. caspase or Annexin V assays).
- The PDF parser expects the standard ProTox-3 result layout. If ProTox changes its page layout, the table conversion may need an update.
- The CB-Dock2 file must contain the columns `ID, Cavities_volume, center_x, center_y, center_z, size_x, size_y, size_z, score, contact_resi`.

## Credits

ProTox-3 and CB-Dock2 are external tools by their respective authors. Please cite them in any report that uses their results.

## License

Add a license of your choice (for example MIT) before publishing.
