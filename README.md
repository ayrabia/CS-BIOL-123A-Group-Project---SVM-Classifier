# Machine Learning-Based Identification of Small Molecule Inhibitors of Mutant BRAF (V600E)

**CS/BIOL 123A — Bioinformatics | SJSU | Spring 2026**

---

## Project Description

This project applies machine learning to identify novel small molecule candidates capable of inhibiting the BRAF V600E oncogenic mutation — one of the most clinically significant driver mutations in melanoma, colorectal cancer, and thyroid cancer. BRAF V600E constitutively activates the MAPK/ERK signaling pathway, leading to uncontrolled cell proliferation. While FDA-approved inhibitors such as vemurafenib exist, acquired drug resistance motivates the search for new candidates.

We build a **Support Vector Machine classifier (SVC — classification, not SVR)** using scikit-learn, trained on instructor-provided CSV datasets of candidate compounds. Molecular descriptors and Morgan fingerprints are computed with RDKit. The trained model is then applied to screen novel candidate compounds retrieved from PubChem. Each predicted inhibitor is assessed by model metrics (precision, recall, accuracy) and cross-referenced with PharmGKB for pharmacogenomic annotations, known drug interactions, and existing clinical data.

---

## Pipeline Overview

```
Instructor CSV Data
(CID, Class label, ~300 precomputed molecular descriptors)
        |
        v
  Task 1 — Data Collection & Cleaning
  (filter records, standardize active/inactive labels
   based on IC50 thresholds; ChEMBL for background context)
        |
        v
  Task 2 — Feature Engineering
  (RDKit: Morgan fingerprints, MW, LogP, HBD, HBA, etc.)
        |
        v
  Task 3 — Exploratory Data Analysis
  (descriptor distributions, class balance, correlations;
   remove low-variance and highly correlated features)
        |
        v
  Task 4 — SVC Model Development
  (scikit-learn SVC; optional Random Forest comparison;
   cross-validation and hyperparameter tuning)
        |
        v
  Task 5 — Model Evaluation
  (Accuracy, Precision, Recall — required;
   AUC-ROC and MCC for imbalanced-class context)
        |
        v
  Task 6 — PubChem Screening
  (retrieve BRAF V600E-relevant compounds via PubChem API;
   compute RDKit descriptors; apply trained SVC to predict
   novel potential inhibitors)
        |
        v
  Task 7 — PharmGKB Integration
  (cross-reference each predicted active compound with
   PharmGKB; assess pharmacogenomic relevance, drug
   interactions, and existing annotations per molecule)
        |
        v
  Task 8 — Interpretation & Reporting
  (identify top predictive features; summarize findings;
   compile final report and code documentation)
```

---

## Team

| Name | Major | Email |
|------|-------|-------|
| Ayman Rabia | Data Science | ayman.rabia@gmail.com |
| Jasjit Deol | Data Science | jasjit.deol@sjsu.edu |
| Chandrahas Sai Vishnumolakala | Data Science | chandrahassai.vishnumolakala@sjsu.edu |
| *(Team Member 4 — TBD)* | — | — |

---

## Repository Structure

```
CS123A_Project/
├── data/
│   ├── raw/                        # Instructor-provided CSV datasets
│   └── processed/                  # Cleaned and feature-engineered data
├── notebooks/
│   ├── 01_data_cleaning.ipynb      # Task 1: load, clean, label compounds
│   ├── 02_feature_engineering.ipynb # Task 2: RDKit descriptors & fingerprints
│   ├── 03_eda.ipynb                # Task 3: distributions, balance, correlations
│   ├── 04_model_training.ipynb     # Task 4: SVC training and tuning
│   ├── 05_model_evaluation.ipynb   # Task 5: metrics — accuracy, precision, recall, AUC-ROC, MCC
│   ├── 06_pubchem_screening.ipynb  # Task 6: PubChem retrieval + SVC predictions
│   ├── 07_pharmgkb_integration.ipynb # Task 7: PharmGKB cross-reference per compound
│   └── 08_interpretation.ipynb     # Task 8: feature analysis, summary, report
├── src/
│   ├── features.py     # RDKit descriptor/fingerprint computation
│   ├── train.py        # SVC training, grid search, model persistence
│   ├── evaluate.py     # accuracy, precision, recall, AUC-ROC, MCC
│   ├── screen.py       # PubChem compound retrieval + SVC screening
│   └── pharmgkb.py     # PharmGKB annotation lookup per compound
├── results/
│   ├── figures/                    # EDA plots, ROC curves, feature importance
│   └── predictions/                # SVC predictions on PubChem candidates
├── requirements.txt
└── README.md
```

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/ayrabia/CS-BIOL-123A-Group-Project---SVM-Classifier.git
cd CS-BIOL-123A-Group-Project---SVM-Classifier
```

### 2. Create a virtual environment

```bash
python3 -m venv venv
source venv/bin/activate   # Windows: venv\Scripts\activate
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

> **Note on RDKit:** RDKit is best installed via conda rather than pip, as the pip distribution can be unstable on some platforms.
> If you encounter issues with `pip install rdkit`, use conda instead:
>
> ```bash
> conda create -n braf-ml python=3.10
> conda activate braf-ml
> conda install -c conda-forge rdkit
> pip install -r requirements.txt   # installs remaining packages
> ```

---

## Data

The primary training dataset is an instructor-provided CSV (`chemical_compounds.csv`) containing PubChem CIDs, activity class labels (0 = inactive, 1 = active), and ~300 precomputed molecular descriptors (PubChem properties and MOE-style descriptors). Place it in `data/raw/` before running the notebooks.

> Note: the training CSV already contains precomputed descriptors, so RDKit is used primarily for featurizing novel PubChem compounds during the screening step (Task 6).

Candidate compounds for screening are retrieved live from the PubChem REST API via `src/pubchem_api.py`.

---

## How to Run the Pipeline

Run notebooks in order from the `notebooks/` directory:

```bash
jupyter notebook
```

| Step | Notebook | SOW Task |
|------|----------|----------|
| 1 | `01_data_cleaning.ipynb` | Load CSV, filter records, standardize active/inactive labels |
| 2 | `02_feature_engineering.ipynb` | Compute RDKit Morgan fingerprints and molecular descriptors |
| 3 | `03_eda.ipynb` | Visualize distributions, class balance, feature correlations |
| 4 | `04_model_training.ipynb` | Train and tune SVC; optional Random Forest comparison |
| 5 | `05_model_evaluation.ipynb` | Report accuracy, precision, recall, AUC-ROC, MCC |
| 6 | `06_pubchem_screening.ipynb` | Retrieve PubChem compounds, compute descriptors, run SVC |
| 7 | `07_pharmgkb_integration.ipynb` | Cross-reference each predicted compound with PharmGKB |
| 8 | `08_interpretation.ipynb` | Feature analysis, findings summary, report compilation |

---

## Classifier Design

We use **scikit-learn's `SVC`** (Support Vector Classifier — the classification variant, not `SVR`) to build a binary small molecule activity prediction model. The model distinguishes active BRAF V600E inhibitors from inactive compounds.

Prediction quality for each novel compound is assessed by:
1. **Model metrics** — precision, recall, and accuracy of the trained SVC
2. **PharmGKB lookup** — if the compound exists in PharmGKB, existing pharmacogenomic annotations, drug interactions, and clinical data are used to further evaluate its relevance

---

## Resources

**Databases:**
- Instructor-provided CSV — primary training dataset (compound CIDs + bioactivity labels)
- PubChem — screening library for novel candidate compounds
- PharmGKB — pharmacogenomic annotations for predicted actives
- UniProt — BRAF protein sequence and functional annotation
- PDB — 3D structure of BRAF (contextual reference)
- ChEMBL — background reference for BRAF bioactivity data

**Software:**
- Python 3, Jupyter Notebook
- RDKit — molecular descriptor computation and fingerprint generation
- scikit-learn — SVC training, cross-validation, evaluation
- pandas, numpy — data manipulation
- matplotlib, seaborn — visualization
- requests — PubChem and PharmGKB API access

---

## Evaluation Metrics

| Metric | Notes |
|--------|-------|
| Accuracy | Required |
| Precision | Required |
| Recall | Required |
| AUC-ROC | Additional; important for imbalanced datasets |
| MCC | Additional; most informative metric under class imbalance |

Class imbalance (many more inactive than active compounds) is expected and addressed via appropriate metric selection and potential class weighting in the SVC.

---

## Schedule

| Weeks | Tasks |
|-------|-------|
| 1–2 | Project setup, literature review, data loading and cleaning (Task 1) |
| 3 | Feature engineering with RDKit, EDA (Tasks 2–3) |
| 4 | SVC model development and initial training (Task 4) |
| 5 | Model tuning, evaluation, and comparison (Task 5) |
| 5–6 | PubChem compound retrieval and SVC-based screening (Task 6) |
| 6 | PharmGKB cross-referencing and pharmacogenomic analysis (Task 7) |
| 7 | Results interpretation, report writing, code cleanup (Task 8) |
| Final | Finalize report, submission package, peer review |

---

## References

- BRAF V600E and papillary thyroid carcinoma outcomes — *J. Clinical Endocrinology & Metabolism*
- BRAF V600E clinicopathological features meta-analysis — *J. Clinical Endocrinology & Metabolism*
- 3D-QSAR and SVM prediction of BRAF-V600E and HIV Integrase inhibitors — instructor-provided publication
- Targeting BRAF V600E in Multiple Myeloma — *Cancer Discovery*
- PubChem Substance and Compound databases — *Nucleic Acids Research*
- PharmGKB — [pharmgkb.org](https://www.pharmgkb.org)
- RDKit — [rdkit.org](https://www.rdkit.org)
