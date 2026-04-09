"""
train.py — SVC model training and persistence.
"""

import joblib
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline


def build_pipeline(C: float = 1.0, kernel: str = "rbf", class_weight: str = "balanced") -> Pipeline:
    """Return a StandardScaler + SVC pipeline."""
    return Pipeline([
        ("scaler", StandardScaler()),
        ("svc", SVC(C=C, kernel=kernel, class_weight=class_weight, probability=True)),
    ])


def tune_and_train(X_train, y_train, cv: int = 5) -> Pipeline:
    """Grid search over C and kernel, return the best fitted pipeline."""
    param_grid = {
        "svc__C":      [0.1, 1, 10, 100],
        "svc__kernel": ["rbf", "linear"],
    }
    base = build_pipeline()
    cv_strategy = StratifiedKFold(n_splits=cv, shuffle=True, random_state=42)
    search = GridSearchCV(base, param_grid, cv=cv_strategy, scoring="roc_auc", n_jobs=-1)
    search.fit(X_train, y_train)
    print(f"Best params: {search.best_params_}  |  CV AUC-ROC: {search.best_score_:.4f}")
    return search.best_estimator_


def save_model(model: Pipeline, path: str) -> None:
    joblib.dump(model, path)
    print(f"Model saved to {path}")


def load_model(path: str) -> Pipeline:
    return joblib.load(path)
