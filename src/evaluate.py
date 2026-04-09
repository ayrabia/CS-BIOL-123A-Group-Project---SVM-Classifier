"""
evaluate.py — Model evaluation metrics (Task 5).

Required: accuracy, precision, recall.
Additional: AUC-ROC, Matthews Correlation Coefficient (MCC).
"""

from sklearn.metrics import (
    accuracy_score,
    precision_score,
    recall_score,
    roc_auc_score,
    matthews_corrcoef,
    classification_report,
    RocCurveDisplay,
)
import matplotlib.pyplot as plt


def report_metrics(y_true, y_pred, y_prob=None) -> dict:
    """Print and return all required and additional metrics."""
    metrics = {
        "accuracy":  accuracy_score(y_true, y_pred),
        "precision": precision_score(y_true, y_pred, zero_division=0),
        "recall":    recall_score(y_true, y_pred, zero_division=0),
        "mcc":       matthews_corrcoef(y_true, y_pred),
    }
    if y_prob is not None:
        metrics["auc_roc"] = roc_auc_score(y_true, y_prob)

    print("\n--- Evaluation Metrics ---")
    for name, val in metrics.items():
        print(f"  {name.upper():12s}: {val:.4f}")
    print("\nClassification Report:")
    print(classification_report(y_true, y_pred, target_names=["Inactive", "Active"]))
    return metrics


def plot_roc(model, X_test, y_test, save_path: str | None = None) -> None:
    """Plot and optionally save an ROC curve."""
    fig, ax = plt.subplots()
    RocCurveDisplay.from_estimator(model, X_test, y_test, ax=ax)
    ax.set_title("SVC ROC Curve — BRAF V600E Inhibitor Classification")
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"ROC curve saved to {save_path}")
    plt.show()
