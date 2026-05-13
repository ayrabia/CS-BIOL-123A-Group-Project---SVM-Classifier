CS/BIOL 123A Project - Spring 2026
BRAF V600E Inhibitor Identification using SVM

Team: Ayman Rabia, Jasjit Deol, Chandrahas Sai Vishnumolakala

------------------------------------------------------------
SETUP
------------------------------------------------------------

Python 3.10+ required.

    pip install -r requirements.txt

If rdkit fails to install via pip, use conda:

    conda create -n braf-ml python=3.10
    conda activate braf-ml
    conda install -c conda-forge rdkit
    pip install -r requirements.txt

------------------------------------------------------------
DATA
------------------------------------------------------------

Place the instructor-provided CSV in data/raw/ before running:

    data/raw/chemical_compounds.csv

All processed files (scaled arrays, model pkl files, predictions)
are generated automatically when you run the notebooks in order.

------------------------------------------------------------
RUNNING THE PIPELINE
------------------------------------------------------------

Launch Jupyter from the project root (not from inside notebooks/):

    jupyter notebook

Run notebooks in this order:

    01_data_cleaning.ipynb       load CSV, remove non-numeric/zero-variance
                                 columns, impute missing values → 311 features
    02_feature_engineering.ipynb scale features with StandardScaler, save arrays
    03_eda.ipynb                 correlation pruning 311→206 features,
                                 Mann-Whitney U feature ranking
    04_model_training.ipynb      GridSearchCV (20 configs, 5-fold CV),
                                 best: C=10, gamma=0.001, RBF kernel
                                 saves svc_model.pkl and scaler.pkl
    05_model_evaluation.ipynb    accuracy/precision/recall/F1/AUC-ROC/MCC
                                 all = 1.0 on held-out test set (n=49)
    06_pubchem_screening.ipynb   similarity search → 591 candidates,
                                 RDKit SVC predicts 533 active (needs internet)
    07_pharmgkb_integration.ipynb cross-references 533 predicted actives
                                 with PharmGKB, finds dabrafenib (needs internet)
    08_interpretation.ipynb      top features, novel hit ranking, summary

------------------------------------------------------------
PROJECT FILES
------------------------------------------------------------

src/
    features.py   RDKit Morgan fingerprint + descriptor computation
    train.py      SVC GridSearchCV training pipeline
    evaluate.py   metric calculations and ROC plotting
    screen.py     PubChem similarity search + screening SVC pipeline
    pharmgkb.py   PharmGKB synonym-mediated lookup

data/raw/          place chemical_compounds.csv here
data/processed/    generated: compounds_clean.csv, X_scaled.npy,
                   y.npy, cid.npy, feature_names.npy,
                   svc_model.pkl, scaler.pkl
results/figures/   all plots and table images
results/predictions/  screening_predictions.csv, predicted_actives.csv,
                      pharmgkb_annotated.csv, pharmgkb_matches_detail.csv

------------------------------------------------------------
KEY RESULTS
------------------------------------------------------------

Dataset:       243 compounds (171 active, 72 inactive), 357 raw columns
After cleaning: 311 features; after EDA: 206 features
Model:         SVC, C=10, gamma=0.001, RBF, class_weight=balanced
Test metrics:  accuracy=1.0, precision=1.0, recall=1.0,
               F1=1.0, AUC-ROC=1.0, MCC=1.0
10-fold CV F1: 0.9957
Screening:     591 PubChem candidates → 533 predicted active
PharmGKB:      1 match — dabrafenib (CID 44462760, FDA-approved BRAF inhibitor)

------------------------------------------------------------
NOTES
------------------------------------------------------------

- Run notebooks from the project root so relative paths (../data/) work
- Notebooks 06 and 07 call live APIs; results may vary slightly by run date
- The screening SVC uses RDKit features (Morgan ECFP4 + 6 descriptors),
  separate from the primary MOE-descriptor-based model
