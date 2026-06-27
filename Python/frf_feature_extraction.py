import pandas as pd
import numpy as np
from scipy.stats import skew, kurtosis

BAND = (0, 400)
N_BANDS = 8


def features_from_spectrum(f, A):
    eps = 1e-12
    feats = {}

    feats["rms"] = np.sqrt(np.mean(A**2))

    peak_i = int(np.argmax(A))
    feats["peak_freq"] = f[peak_i]
    feats["peak_mag"] = A[peak_i]

    centroid = np.sum(f * A) / (np.sum(A) + eps)
    feats["centroid"] = centroid

    feats["spread"] = np.sqrt(
        np.sum(((f - centroid) ** 2) * A) / (np.sum(A) + eps)
    )

    feats["skew"] = skew(A)
    feats["kurt"] = kurtosis(A)

    edges = np.linspace(BAND[0], BAND[1], N_BANDS + 1)

    for b in range(N_BANDS):
        sel = (f >= edges[b]) & (f < edges[b + 1])
        feats[f"band_{b}"] = np.sum(A[sel] ** 2)

    return feats


def build_allowed_runs():
    """
    Creates the exact dataset according to your measurement protocol.

    Healthy:
    H1, H2, H3
    P1-P4
    S1-S4
    = 3 x 4 x 4 = 48 healthy samples

    Abnormal:
    D1 P1-P4 S1-S4 = 16 damaged samples
    H1P2S2NA screw-removal cases = 6 samples
    D1P2S2DA screw-removal cases = 6 samples
    Total abnormal = 28 samples
    """

    healthy_runs = set()

    for specimen in ["H1", "H2", "H3"]:
        for p in range(1, 5):
            for s in range(1, 5):
                healthy_runs.add(f"{specimen}P{p}S{s}")

    abnormal_runs = set()

    # Built-in damaged tube: 4 positions x 4 strikes = 16
    for p in range(1, 5):
        for s in range(1, 5):
            abnormal_runs.add(f"D1P{p}S{s}")

    # Healthy tube with screw removals at P2S2
    healthy_screw_cases = [
        "H1P2S2NA1",
        "H1P2S2NA2",
        "H1P2S2NA3",
        "H1P2S2NA4",
        "H1P2S2NA1NA2",
        "H1P2S2NA3NA4",
    ]

    # Damaged tube with screw removals at P2S2
    damaged_screw_cases = [
        "D1P2S2DA1",
        "D1P2S2DA2",
        "D1P2S2DA3",
        "D1P2S2DA4",
        "D1P2S2DA1DA2",
        "D1P2S2DA3DA4",
    ]

    abnormal_runs.update(healthy_screw_cases)
    abnormal_runs.update(damaged_screw_cases)

    return healthy_runs, abnormal_runs


def assign_label(run_name):
    healthy_runs, abnormal_runs = build_allowed_runs()

    if run_name in healthy_runs:
        return 0

    if run_name in abnormal_runs:
        return 1

    return None


def main():
    df = pd.read_csv("Section1_FRF_only/Section1_FRF_combined_long.csv")

    df["run"] = df["run"].astype(str).str.strip()
    df["sensor"] = df["sensor"].astype(str).str.strip()

    rows = []
    audit_rows = []

    for run_name, run_df in df.groupby("run"):
        label = assign_label(run_name)

        if label is None:
            audit_rows.append({
                "run": run_name,
                "status": "excluded",
                "label": None
            })
            continue

        row = {}

        for sensor_name, sensor_df in run_df.groupby("sensor"):
            sensor_df = sensor_df.sort_values("x_value")

            f = sensor_df["x_value"].values
            A = sensor_df["magnitude"].values

            feats = features_from_spectrum(f, A)

            for key, value in feats.items():
                row[f"{sensor_name}_{key}"] = value

        row["label"] = label
        row["campaign_id"] = run_name

        rows.append(row)

        audit_rows.append({
            "run": run_name,
            "status": "included",
            "label": label
        })

    features = pd.DataFrame(rows)
    audit = pd.DataFrame(audit_rows)

    features.to_csv("features_frf.csv", index=False)
    audit.to_csv("frf_labels_audit.csv", index=False)

    print("Saved features_frf.csv")
    print("Saved frf_labels_audit.csv")
    print(features.shape)
    print(features["label"].value_counts())


if __name__ == "__main__":
    main()