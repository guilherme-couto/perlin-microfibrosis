import os
import pandas as pd
import numpy as np
import joblib
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

# ========== Settings ==========
input_file = "./fibrosis_metrics.csv"
output_dir = "./pca_outputs"
os.makedirs(output_dir, exist_ok=True)
for mode in ["threshold", "composition", "general"]:
    os.makedirs(f"{output_dir}/plots/{mode}", exist_ok=True)
    os.makedirs(f"{output_dir}/data/{mode}", exist_ok=True)

# ========== Load data ==========
df_all = pd.read_csv(input_file)

# ========== Define reference densities ==========
reference_densities = {
    "interstitial": 0.096,
    "compact": 0.472,
    "diffuse": 0.22,
    "patchy": 0.269
}

# ========== Prepare data ==========
feature_cols = [col for col in df_all.columns if any(col.startswith(prefix) for prefix in 
                                                     ["orientation_", "major_axis_", "minor_axis_"])]
print(f"Feature columns used for PCA total({len(feature_cols)}):", feature_cols)

df_all["fibro_type"] = df_all["fibro_typename"]
df_all["mode"] = df_all["generation_mode"]
df_all["density"] = df_all["density"].astype(float)
theta = df_all["theta"]
y_type = df_all["fibro_type"]
y_mode = df_all["mode"]

# ========== Standardize reference features ==========
# Check if scaler already exists
if os.path.exists(f"{output_dir}/data/threshold/scaler.pkl"):
    scaler = joblib.load(f"{output_dir}/data/threshold/scaler.pkl")

else:
    # ========== Select reference data for PCA ==========
    df_ref = pd.concat([
        df_all[
            (df_all["mode"] == "threshold") &
            (df_all["fibro_type"] == ftype) &
            (np.isclose(df_all["density"], ref_density))
        ]
        for ftype, ref_density in reference_densities.items()
    ], ignore_index=True)

    X_ref = df_ref[feature_cols].values
    y_ref_type = df_ref["fibro_type"].values
    theta_ref = df_ref["theta"].values

    print("Creating new scaler for reference data.")
    scaler = StandardScaler()
    X_ref_scaled = scaler.fit_transform(X_ref)

# ========== PCA from reference data ==========
# Check if PCA model already exists
if os.path.exists(f"{output_dir}/data/threshold/pca_2d.pkl"):
    pca_2d = joblib.load(f"{output_dir}/data/threshold/pca_2d.pkl")

else:
    print("Creating new PCA model for reference data.")
    pca_2d = PCA(n_components=2)
    X_ref_pca_2d = pca_2d.fit_transform(X_ref_scaled)

    # ========== Save scaler and PCA model ==========
    print("Saving scaler and PCA model.")
    joblib.dump(scaler, f"{output_dir}/data/threshold/scaler.pkl")
    joblib.dump(pca_2d, f"{output_dir}/data/threshold/pca_2d.pkl")

    # Save explained variance
    explained_variance = pca_2d.explained_variance_ratio_
    np.savetxt(f"{output_dir}/data/threshold/explained_variance_ratio.txt", explained_variance, fmt="%.6f")
    print(f"Explained variance ratio (PC1, PC2): {explained_variance}")

    # Save reference projection
    df_ref_pca_2d = pd.DataFrame(X_ref_pca_2d, columns=["PC1", "PC2"])
    df_ref_pca_2d["fibro_type"] = y_ref_type
    df_ref_pca_2d["theta"] = theta_ref
    df_ref_pca_2d.to_csv(f"{output_dir}/data/threshold/pca_2d_ref.csv", index=False)

# Project all data using reference PCA
X_all_scaled = scaler.transform(df_all[feature_cols].values)
X_all_pca_2d = pca_2d.transform(X_all_scaled)

df_all_pca_2d = pd.DataFrame(X_all_pca_2d, columns=["PC1", "PC2"])
df_all_pca_2d["fibro_type"] = y_type
df_all_pca_2d["theta"] = theta
df_all_pca_2d["mode"] = y_mode
df_all_pca_2d.to_csv(f"{output_dir}/data/general/pca_2d_all.csv", index=False)

# ========== Plots ==========
def plot_by_theta(df_subset, title, save_path):
    plt.figure(figsize=(10, 7))
    scatter = plt.scatter(df_subset["PC1"], df_subset["PC2"], c=df_subset["theta"], cmap="viridis", s=70)
    plt.colorbar(scatter, label="Theta (radians)")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title(title)
    plt.savefig(save_path, dpi=300)
    plt.close()

# Plot PCA colored by theta
plot_by_theta(df_all_pca_2d[df_all_pca_2d["mode"] == "threshold"],
              "PCA (2D) of Threshold Patterns Colored by Theta",
              f"{output_dir}/plots/threshold/pca_2d_theta.png")

plot_by_theta(df_all_pca_2d[df_all_pca_2d["mode"] == "composition"],
              "PCA (2D) of Composition Patterns Colored by Theta",
              f"{output_dir}/plots/composition/pca_2d_theta.png")

plot_by_theta(df_all_pca_2d,
              "PCA (2D) of All Patterns Colored by Theta",
              f"{output_dir}/plots/general/pca_2d_theta.png")

# Existing plot by type and generation mode
plt.figure(figsize=(10, 7))
sns.scatterplot(data=df_all_pca_2d[df_all_pca_2d["mode"] == "threshold"],
                x="PC1", y="PC2", hue="fibro_type", palette="tab10", s=70)
plt.title("PCA (2D) of Threshold Fibrosis Patterns by Type")
plt.savefig(f"{output_dir}/plots/threshold/pca_2d_by_type.png", dpi=300)
plt.close()

plt.figure(figsize=(10, 7))
sns.scatterplot(data=df_all_pca_2d[df_all_pca_2d["mode"] == "composition"],
                x="PC1", y="PC2", hue="fibro_type", palette="tab10", s=70)
plt.title("PCA (2D) of Composition Fibrosis Patterns by Type")
plt.savefig(f"{output_dir}/plots/composition/pca_2d_by_type.png", dpi=300)
plt.close()

plt.figure(figsize=(10, 7))
sns.scatterplot(data=df_all_pca_2d, x="PC1", y="PC2", hue="fibro_type", style="mode", palette="tab10", s=70)
plt.title("PCA (2D) of All Fibrosis Patterns by Type and Generation Mode")
plt.savefig(f"{output_dir}/plots/general/pca_2d.png", dpi=300)
plt.close()

# ========== KMeans Clustering (optional) ==========
kmeans = KMeans(n_clusters=4, random_state=42)
df_all_pca_2d["cluster"] = kmeans.fit_predict(X_all_pca_2d)
df_all_pca_2d.to_csv(f"{output_dir}/data/general/pca_2d_with_clusters.csv", index=False)

plt.figure(figsize=(10, 7))
sns.scatterplot(data=df_all_pca_2d[df_all_pca_2d["mode"] == "threshold"], x="PC1", y="PC2", hue="cluster", palette="Dark2", s=70)
plt.title("PCA (2D) of Threshold Patterns with KMeans Clusters")
plt.savefig(f"{output_dir}/plots/threshold/pca_2d_clusters.png", dpi=300)
plt.close()

plt.figure(figsize=(10, 7))
sns.scatterplot(data=df_all_pca_2d[df_all_pca_2d["mode"] == "composition"], x="PC1", y="PC2", hue="cluster", palette="Dark2", s=70)
plt.title("PCA (2D) of Composition Patterns with KMeans Clusters")
plt.savefig(f"{output_dir}/plots/composition/pca_2d_clusters.png", dpi=300)
plt.close()

plt.figure(figsize=(10, 7))
sns.scatterplot(data=df_all_pca_2d, x="PC1", y="PC2", hue="cluster", style="mode", palette="Dark2", s=70)
plt.title("PCA (2D) of All Patterns with KMeans Clusters")
plt.savefig(f"{output_dir}/plots/general/pca_2d_clusters.png", dpi=300)
plt.close()

print("PCA complete. Outputs saved in:", output_dir)
