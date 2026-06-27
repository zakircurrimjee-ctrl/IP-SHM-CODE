import warnings
import numpy as np
import pandas as pd

from sklearn.exceptions import UndefinedMetricWarning
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, f1_score
from sklearn.pipeline import Pipeline

warnings.filterwarnings("ignore", category=UndefinedMetricWarning)

# Load the FRF feature dataset
df = pd.read_csv("features_frf.csv")

# Remove rows with missing values, just in case
df = df.dropna()

feature_cols = [
    "Acc18_band_1",
    "Acc18_peak_freq",
]
print("Using selected features:")
for col in feature_cols:
    print("-", col)

X = df[feature_cols].values
y = df["label"].values

print("Dataset shape:", df.shape)
print("Labels:")
print(df["label"].value_counts())

# 70/30 train-test split, stratified by label
X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.30,
    random_state=42,
    stratify=y
)

def report_model(name, model):
    pred = model.predict(X_test)

    print(f"\n=== {name} ===")
    print("Confusion matrix [[TN FP] [FN TP]]:")
    print(confusion_matrix(y_test, pred))

    print(f"Accuracy : {accuracy_score(y_test, pred):.3f}")
    print(f"Precision: {precision_score(y_test, pred, zero_division=0):.3f}")
    print(f"Recall   : {recall_score(y_test, pred, zero_division=0):.3f}")
    print(f"F1 score : {f1_score(y_test, pred, zero_division=0):.3f}")

# -------------------------
# Random Forest
# -------------------------
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

print("\nBest Random Forest parameters:")
print(rf_search.best_params_)
print("Best Random Forest cross-validation F1:", round(rf_search.best_score_, 3))

report_model("Random Forest", rf_search)

# Feature importance
best_rf = rf_search.best_estimator_.named_steps["clf"]
importance = pd.Series(
    best_rf.feature_importances_,
    index=feature_cols
).sort_values(ascending=False)

print("\nRandom Forest feature importance:")
print(importance.round(3).to_string())

# -------------------------
# SVM
# -------------------------
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

print("\nBest SVM parameters:")
print(svm_search.best_params_)
print("Best SVM cross-validation F1:", round(svm_search.best_score_, 3))
report_model("SVM RBF", svm_search)
def report_model(name, model):
    train_pred = model.predict(X_train)
    test_pred = model.predict(X_test)

    print(f"\n=== {name} ===")

    print("\nTraining set:")
    print("Confusion matrix [[TN FP] [FN TP]]:")
    print(confusion_matrix(y_train, train_pred))
    print(f"Accuracy : {accuracy_score(y_train, train_pred):.3f}")
    print(f"Precision: {precision_score(y_train, train_pred, zero_division=0):.3f}")
    print(f"Recall   : {recall_score(y_train, train_pred, zero_division=0):.3f}")
    print(f"F1 score : {f1_score(y_train, train_pred, zero_division=0):.3f}")

    print("\nTest set:")
    print("Confusion matrix [[TN FP] [FN TP]]:")
    print(confusion_matrix(y_test, test_pred))
    print(f"Accuracy : {accuracy_score(y_test, test_pred):.3f}")
    print(f"Precision: {precision_score(y_test, test_pred, zero_division=0):.3f}")
    print(f"Recall   : {recall_score(y_test, test_pred, zero_division=0):.3f}")
    print(f"F1 score : {f1_score(y_test, test_pred, zero_division=0):.3f}")