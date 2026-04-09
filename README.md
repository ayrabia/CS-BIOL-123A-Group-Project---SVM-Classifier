# Machine Learning-Based Identification of Small Molecule Inhibitors of Mutant BRAF (V600E)

**CS/BIOL 123A — Bioinformatics | SJSU | Spring 2026**

---

## Project Description

This project applies machine learning to identify novel small molecule candidates capable of inhibiting the BRAF V600E oncogenic mutation — one of the most clinically significant driver mutations in melanoma, colorectal cancer, and thyroid cancer. BRAF V600E constitutively activates the MAPK/ERK signaling pathway, leading to uncontrolled cell proliferation. While FDA-approved inhibitors such as vemurafenib exist, acquired drug resistance motivates the search for new candidates.

We build a Support Vector Machine (SVC) classifier trained on instructor-provided compound datasets with precomputed molecular descriptors. The trained model is then used to screen candidate compounds retrieved from PubChem, and top predicted inhibitors are cross-referenced with PharmGKB for pharmacogenomic context.

---

## Pipeline Overview

```
Instructor CSV Data
(CID, Class, ~300 molecular descriptors)
        |
        v
  Data Cleaning &
  Label Standardization
  (active / inactive based on IC50 threshold)
        |
        v
  Feature Engineering
  (RDKit: Morgan fingerprints, MW, LogP, HBD, HBA)
        |
        v
  Exploratory Data Analysis
  (descriptor distributions, class balance, correlations)
        |
        v
  SVC Classifier Training
  (scikit-learn; optional: Random Forest comparison)
        |
        v
  Model Evaluation
  (Accuracy, Precision, Recall, AUC-ROC, MCC)
        |
        v
  PubChem Screening
  (retrieve BRAF V600E-relevant compounds via API,
   compute descriptors with RDKit, run SVC predictions)
        |
        v
  PharmGKB Validation
  (cross-reference predicted actives for known drug
   interactions & pharmacogenomic annotations)
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
│   ├── raw/                  # Instructor-provided CSV datasets
│   └── processed/            # Cleaned and feature-engineered data
├── notebooks/
│   ├── 01_data_cleaning.ipynb
│   ├── 02_feature_engineering.ipynb
│   ├── 03_eda.ipynb
│   ├── 04_model_training.ipynb
│   ├── 05_pubchem_screening.ipynb
│   └── 06_pharmgkb_validation.ipynb
├── src/
│   ├── features.py           # RDKit descriptor/fingerprint computation
│   ├── model.py              # SVC training and evaluation utilities
│   ├── pubchem_api.py        # PubChem compound retrieval
│   └── pharmgkb_api.py       # PharmGKB cross-referencing
├── results/
│   ├── figures/              # EDA plots, ROC curves, feature importance
│   └── predictions/          # SVC output on PubChem candidates
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

Candidate compounds for screening are retrieved live from the PubChem REST API using `src/pubchem_api.py`.

---

## How to Run the Pipeline

Run notebooks in order from the `notebooks/` directory:

```bash
jupyter notebook
```

| Step | Notebook | Description |
|------|----------|-------------|
| 1 | `01_data_cleaning.ipynb` | Load CSV, clean records, standardize activity labels |
| 2 | `02_feature_engineering.ipynb` | Compute RDKit Morgan fingerprints and descriptors |
| 3 | `03_eda.ipynb` | Visualize distributions, class balance, feature correlations |
| 4 | `04_model_training.ipynb` | Train and tune SVC; optional Random Forest comparison |
| 5 | `05_pubchem_screening.ipynb` | Retrieve PubChem compounds, apply trained SVC |
| 6 | `06_pharmgkb_validation.ipynb` | Cross-reference predictions with PharmGKB |

---

## Key Dependencies

| Package | Purpose |
|---------|---------|
| `scikit-learn` | SVC classifier, cross-validation, evaluation metrics |
| `rdkit` | Morgan fingerprints, molecular descriptors |
| `pandas` / `numpy` | Data manipulation |
| `matplotlib` / `seaborn` | Visualization |
| `requests` | PubChem and PharmGKB API calls |

---

## Evaluation Metrics

Model performance is reported using:
- Accuracy, Precision, Recall
- AUC-ROC
- Matthews Correlation Coefficient (MCC) — particularly important given expected class imbalance between active and inactive compounds

---

## References

- BRAF V600E and papillary thyroid carcinoma outcomes — *J. Clinical Endocrinology & Metabolism*
- 3D-QSAR and SVM prediction of BRAF-V600E inhibitors — instructor-provided publication
- PubChem Substance and Compound databases — *Nucleic Acids Research*
- PharmGKB — pharmacogenomic annotations ([pharmgkb.org](https://www.pharmgkb.org))
- RDKit documentation — [rdkit.org](https://www.rdkit.org)
