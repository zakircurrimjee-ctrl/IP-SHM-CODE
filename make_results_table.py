import pandas as pd

# ============================================================
# Results from compare_feature_counts.py
# ============================================================

results = [
    # Features, Model, Accuracy, Precision, Recall, F1
    [1, "Random Forest", 0.522, 0.333, 0.375, 0.353],
    [1, "SVM RBF",       0.696, 0.667, 0.250, 0.364],

    [2, "Random Forest", 0.826, 0.833, 0.625, 0.714],
    [2, "SVM RBF",       0.957, 0.889, 1.000, 0.941],

    [3, "Random Forest", 0.826, 0.833, 0.625, 0.714],
    [3, "SVM RBF",       0.870, 0.857, 0.750, 0.800],

    [4, "Random Forest", 0.870, 0.857, 0.750, 0.800],
    [4, "SVM RBF",       0.739, 0.667, 0.500, 0.571],

    [5, "Random Forest", 0.870, 0.857, 0.750, 0.800],
    [5, "SVM RBF",       0.652, 0.500, 0.125, 0.200],

    [8, "Random Forest", 0.739, 0.750, 0.375, 0.500],
    [8, "SVM RBF",       0.739, 0.750, 0.375, 0.500],

    [10, "Random Forest", 0.739, 0.750, 0.375, 0.500],
    [10, "SVM RBF",       0.739, 0.750, 0.375, 0.500],

    [15, "Random Forest", 0.739, 0.750, 0.375, 0.500],
    [15, "SVM RBF",       0.826, 1.000, 0.500, 0.667],

    [30, "Random Forest", 0.739, 1.000, 0.250, 0.400],
    [30, "SVM RBF",       0.870, 1.000, 0.625, 0.769],
]

df = pd.DataFrame(
    results,
    columns=["Features used", "Model", "Accuracy", "Precision", "Recall", "F1-score"]
)

# ============================================================
# Learning curve interpretation from your plots
# These are approximate values read from the final point
# at about 60 training samples.
# ============================================================

learning_curve_results = [
    [2, "Random Forest", "~0.83", "Good; validation F1 increases with more data"],
    [4, "Random Forest", "~0.80", "Good; validation F1 improves at larger sample size"],
    [5, "Random Forest", "~0.88", "Best learning-curve behaviour among tested models"],

    [2, "SVM RBF", "~0.33", "Unstable; high test result not supported by learning curve"],
    [4, "SVM RBF", "~0.33", "Weak validation performance"],
    [5, "SVM RBF", "~0.57", "Improves with more data, but still below Random Forest"],
]

lc = pd.DataFrame(
    learning_curve_results,
    columns=[
        "Features used",
        "Model",
        "Approx. final validation F1 from learning curve",
        "Learning-curve interpretation",
    ]
)

df = df.merge(
    lc,
    on=["Features used", "Model"],
    how="left"
)

# ============================================================
# Improvements compared with each model's own 30-feature baseline
# ============================================================

baseline = (
    df[df["Features used"] == 30]
    .set_index("Model")[["Accuracy", "Recall", "F1-score"]]
)

def improvement(row, metric):
    model = row["Model"]
    return row[metric] - baseline.loc[model, metric]

df["F1 improvement vs 30 features"] = df.apply(
    lambda row: improvement(row, "F1-score"), axis=1
)

df["Recall improvement vs 30 features"] = df.apply(
    lambda row: improvement(row, "Recall"), axis=1
)

# ============================================================
# Add clear interpretation
# ============================================================

def add_interpretation(row):
    model = row["Model"]
    k = row["Features used"]

    if model == "Random Forest" and k == 5:
        return (
            "Selected final reduced-feature model: strong test result and best learning-curve behaviour"
        )

    if model == "Random Forest" and k == 4:
        return (
            "Best RF test result, tied with 5 features; learning curve slightly lower than 5-feature RF"
        )

    if model == "Random Forest" and k == 2:
        return (
            "Improves RF compared with 30 features, but not as strong as 4/5 features"
        )

    if model == "SVM RBF" and k == 2:
        return (
            "Best single train/test result, but learning curve suggests unstable generalisation"
        )

    if model == "SVM RBF" and k == 30:
        return (
            "Strong full-feature SVM baseline, but more complex and lower recall than 2-feature single-split result"
        )

    if k == 30:
        return "Full-feature baseline"

    if row["F1 improvement vs 30 features"] > 0:
        return "Improved compared with 30-feature baseline"

    if row["F1 improvement vs 30 features"] == 0:
        return "Same F1 as 30-feature baseline"

    return "Lower than 30-feature baseline"

df["Final interpretation"] = df.apply(add_interpretation, axis=1)

# Round numerical values
numeric_cols = [
    "Accuracy",
    "Precision",
    "Recall",
    "F1-score",
    "F1 improvement vs 30 features",
    "Recall improvement vs 30 features",
]

df[numeric_cols] = df[numeric_cols].round(3)

# ============================================================
# Full sensitivity table
# ============================================================

full_table = df[
    [
        "Features used",
        "Model",
        "Accuracy",
        "Precision",
        "Recall",
        "F1-score",
        "F1 improvement vs 30 features",
        "Recall improvement vs 30 features",
        "Approx. final validation F1 from learning curve",
        "Learning-curve interpretation",
        "Final interpretation",
    ]
].sort_values(["Model", "Features used"])

# ============================================================
# Thesis-ready summary table
# ============================================================

summary_rows = [
    {
        "Purpose": "RF full-feature baseline",
        "Model": "Random Forest",
        "Features used": 30,
    },
    {
        "Purpose": "RF reduced-feature option",
        "Model": "Random Forest",
        "Features used": 2,
    },
    {
        "Purpose": "RF best test result",
        "Model": "Random Forest",
        "Features used": 4,
    },
    {
        "Purpose": "RF selected final model",
        "Model": "Random Forest",
        "Features used": 5,
    },
    {
        "Purpose": "SVM full-feature baseline",
        "Model": "SVM RBF",
        "Features used": 30,
    },
    {
        "Purpose": "SVM best single-split result",
        "Model": "SVM RBF",
        "Features used": 2,
    },
]

summary = pd.DataFrame(summary_rows)

summary = summary.merge(
    df,
    on=["Model", "Features used"],
    how="left"
)

summary = summary[
    [
        "Purpose",
        "Model",
        "Features used",
        "Accuracy",
        "Precision",
        "Recall",
        "F1-score",
        "F1 improvement vs 30 features",
        "Recall improvement vs 30 features",
        "Approx. final validation F1 from learning curve",
        "Learning-curve interpretation",
        "Final interpretation",
    ]
]

# Save tables
full_table.to_csv("feature_count_sensitivity_table.csv", index=False)
summary.to_csv("New_model_comparison_table.csv", index=False)

print("\n=== Full feature-count sensitivity table ===")
print(full_table.to_string(index=False))

print("\n=== Thesis-ready model comparison table ===")
print(summary.to_string(index=False))

print("\nSaved:")
print("- feature_count_sensitivity_table.csv")
print("- New_model_comparison_table.csv")