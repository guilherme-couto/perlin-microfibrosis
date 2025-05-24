import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from mpl_toolkits.mplot3d import Axes3D

# ========== Settings ==========
input_dir = "./"  # folder where the CSV files are located
output_dir = "./pca_outputs"
os.makedirs(output_dir, exist_ok=True)
os.makedirs(f"{output_dir}/plots", exist_ok=True)
os.makedirs(f"{output_dir}/data", exist_ok=True)

# ========== Load and combine datasets ==========
filenames = {
    "compact": "metrics_compact.csv",
    "diffuse": "metrics_diffuse.csv",
    "interstitial": "metrics_interstitial.csv",
    "patchy": "metrics_patchy.csv"
}

data_frames = []

for fibro_type, filename in filenames.items():
    path = os.path.join(input_dir, filename)
    df = pd.read_csv(path)
    df["fibro_type"] = fibro_type
    data_frames.append(df)

df_all = pd.concat(data_frames, ignore_index=True)
df_all.to_csv(f"{output_dir}/data/merged_metrics.csv", index=False)

# ========== Prepare data ==========
feature_cols = df_all.columns[2:-1]  # exclude 'fibro_typename', 'seed', 'fibro_type'
if "fibro_typename" in feature_cols:
    feature_cols = feature_cols.drop("fibro_typename")
elif "seed" in feature_cols:
    feature_cols = feature_cols.drop("seed")
elif "fibro_type" in feature_cols:
    feature_cols = feature_cols.drop("fibro_type")

# Ensure feature_cols is a list
feature_cols = list(feature_cols)
print("Feature columns used for PCA:", feature_cols)

X = df_all[feature_cols].values
y = df_all["fibro_type"].values

# ========== Standardize features ==========
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# ========== PCA: 2D ==========
pca_2d = PCA(n_components=2)
X_pca_2d = pca_2d.fit_transform(X_scaled)

df_pca_2d = pd.DataFrame(X_pca_2d, columns=["PC1", "PC2"])
df_pca_2d["fibro_type"] = y
df_pca_2d.to_csv(f"{output_dir}/data/pca_2d.csv", index=False)

# Plot PCA 2D
plt.figure(figsize=(10, 7))
sns.scatterplot(data=df_pca_2d, x="PC1", y="PC2", hue="fibro_type", palette="Set2", s=70)
plt.title("PCA (2D) of Fibrosis Metrics")
plt.savefig(f"{output_dir}/plots/pca_2d.png", dpi=300)
plt.close()

# ========== PCA: 3D ==========
pca_3d = PCA(n_components=3)
X_pca_3d = pca_3d.fit_transform(X_scaled)

df_pca_3d = pd.DataFrame(X_pca_3d, columns=["PC1", "PC2", "PC3"])
df_pca_3d["fibro_type"] = y
df_pca_3d.to_csv(f"{output_dir}/data/pca_3d.csv", index=False)

# Plot PCA 3D
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection="3d")
for label in np.unique(y):
    subset = df_pca_3d[df_pca_3d["fibro_type"] == label]
    ax.scatter(subset["PC1"], subset["PC2"], subset["PC3"], label=label, s=60)
ax.set_title("PCA (3D) of Fibrosis Metrics")
ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
ax.set_zlabel("PC3")
ax.legend()
plt.savefig(f"{output_dir}/plots/pca_3d.png", dpi=300)
plt.close()

# ========== KMeans Clustering on PCA 2D ==========
kmeans = KMeans(n_clusters=4, random_state=42)
df_pca_2d["cluster"] = kmeans.fit_predict(X_pca_2d)
df_pca_2d.to_csv(f"{output_dir}/data/pca_2d_with_clusters.csv", index=False)

# Plot clusters
plt.figure(figsize=(10, 7))
sns.scatterplot(data=df_pca_2d, x="PC1", y="PC2", hue="cluster", style="fibro_type", palette="Set1", s=70)
plt.title("PCA (2D) with KMeans Clustering")
plt.savefig(f"{output_dir}/plots/pca_2d_clusters.png", dpi=300)
plt.close()

print("Analysis complete. Outputs saved to:", output_dir)
