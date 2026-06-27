import pandas as pd
import numpy as np

from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score

df = pd.read_csv("features_frf.csv")
df = df.dropna()

# Feature ranking from your reduced Random Forest result
ranked_features = [
    "Acc18_band_1",
    "Acc18_peak_freq",
    "Acc23_band_1",
    "Acc23_peak_freq",
    "Acc18_kurt",
    "Acc23_band_0",
    "Acc18_band_0",
    "Acc23_kurt",
    "Acc18_band_3",
    "Acc23_band_3",
    "Acc18_skew",
    "Acc18_band_5",
    "Acc23_centroid",
    "Acc18_band_6",
    "Acc18_centroid",
    "Acc23_band_6",
    "Acc23_skew",
    "Acc18_band_2",
    "Acc18_spread",
    "Acc23_band_7",
    "Acc23_band_2",
    "Acc18_peak_mag",
    "Acc18_band_7",
    "Acc23_spread",
    "Acc23_band_5",
    "Acc18_band_4",
    "Acc23_rms",
    "Acc23_band_4",
    "Acc18_rms",
    "Acc23_peak_mag",
]

feature_counts = [1, 2, 3, 4, 5, 8, 10, 15, 30]

results = []

for k in feature_counts:
    feature_cols = ranked_features[:k]

    X = df[feature_cols].values
    y = df["label"].values

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=0.30,
        random_state=42,
        stratify=y
    )

    # Random Forest
    rf_grid = {
        "clf__n_estimators": [100, 300, 500],
        "clf__max_depth": [None, 5, 10],
        "clf__min_samples_leaf": [1, 2, 4],
    }

    rf_pipeline = Pipeline([
        ("scaler", StandardScaler()),
        ("clf", RandomForestClassifier(random_state=42))
    ])

    rf_search = GridSearchCV(
        rf_pipeline,
        rf_grid,
        cv=5,
        scoring="f1"
    )

    rf_search.fit(X_train, y_train)
    rf_pred = rf_search.predict(X_test)

    results.append({
        "Features": k,
        "Model": "Random Forest",
        "Accuracy": accuracy_score(y_test, rf_pred),
        "Precision": precision_score(y_test, rf_pred, zero_division=0),
        "Recall": recall_score(y_test, rf_pred, zero_division=0),
        "F1": f1_score(y_test, rf_pred, zero_division=0),
        "Best parameters": rf_search.best_params_,
    })

    # SVM
    svm_grid = {
        "clf__C": [0.1, 1, 10, 100],
        "clf__gamma": [0.001, 0.01, 0.1, 1],
    }

    svm_pipeline = Pipeline([
        ("scaler", StandardScaler()),
        ("clf", SVC(kernel="rbf"))
    ])

    svm_search = GridSearchCV(
        svm_pipeline,
        svm_grid,
        cv=5,
        scoring="f1"
    )

    svm_search.fit(X_train, y_train)
    svm_pred = svm_search.predict(X_test)

    results.append({
        "Features": k,
        "Model": "SVM RBF",
        "Accuracy": accuracy_score(y_test, svm_pred),
        "Precision": precision_score(y_test, svm_pred, zero_division=0),
        "Recall": recall_score(y_test, svm_pred, zero_division=0),
        "F1": f1_score(y_test, svm_pred, zero_division=0),
        "Best parameters": svm_search.best_params_,
    })

results_df = pd.DataFrame(results)

print(results_df[["Features", "Model", "Accuracy", "Precision", "Recall", "F1"]].round(3).to_string(index=False))

results_df.to_csv("feature_count_comparison.csv", index=False)
print("\nSaved feature_count_comparison.csv")