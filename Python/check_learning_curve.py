import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.model_selection import learning_curve, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier

# Load dataset
df = pd.read_csv("features_frf.csv")
df = df.dropna()

feature_cols = [
    "Acc18_band_1",
    "Acc18_peak_freq",
    "Acc23_band_1",
    "Acc23_peak_freq",
]
X = df[feature_cols].values
y = df["label"].values

print("Dataset shape:", df.shape)
print("Labels:")
print(df["label"].value_counts())

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

def plot_learning_curve(model, model_name, output_name):
    train_sizes, train_scores, val_scores = learning_curve(
        model,
        X,
        y,
        cv=cv,
        scoring="f1",
        train_sizes=[24, 36, 48, 60],
        shuffle=True,
        random_state=42,
        n_jobs=-1
    )

    train_mean = np.mean(train_scores, axis=1)
    train_std = np.std(train_scores, axis=1)

    val_mean = np.mean(val_scores, axis=1)
    val_std = np.std(val_scores, axis=1)

    print(f"\n=== {model_name} Learning Curve ===")
    for size, tr, va in zip(train_sizes, train_mean, val_mean):
        print(f"Training size: {size:>3} | Train F1: {tr:.3f} | Validation F1: {va:.3f}")

    plt.figure()
    plt.plot(train_sizes, train_mean, marker="o", label="Training F1")
    plt.plot(train_sizes, val_mean, marker="o", label="Validation F1")

    plt.fill_between(
        train_sizes,
        train_mean - train_std,
        train_mean + train_std,
        alpha=0.2
    )

    plt.fill_between(
        train_sizes,
        val_mean - val_std,
        val_mean + val_std,
        alpha=0.2
    )

    plt.xlabel("Number of training samples")
    plt.ylabel("F1 score")
    plt.title(f"Learning Curve: {model_name}")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_name, dpi=300)
    plt.close()



# Your best SVM model from the previous result
svm_model = Pipeline([
    ("scaler", StandardScaler()),
    ("clf", SVC(kernel="rbf", C=10, gamma=0.1))
])

# Random Forest model
rf_model = Pipeline([
    ("scaler", StandardScaler()),
    ("clf", RandomForestClassifier(
        n_estimators=500,
        max_depth=None,
        min_samples_leaf=2,
        random_state=42
    ))
])

plot_learning_curve(svm_model, "SVM RBF - 4 Features", "svm_4features_learning_curve.png")
plot_learning_curve(rf_model, "Random Forest - 4 Features", "rf_4features_learning_curve.png")