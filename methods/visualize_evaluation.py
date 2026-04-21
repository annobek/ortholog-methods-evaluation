#import pandas as pd
import matplotlib.pyplot as plt
#import seaborn 
import numpy as np

# Visualize metrics as heatmaps

# Prepare metrics for visualisation

# Vizualize heatmap of methods overview


'''
# for pure cactus metrics
methods_cactus = ["Cactus", 
            "Neighborhood", 
            "Tree-Majority",
            "Tree-Standard",
            "Tree-Whitelist",
            "PTP",
            "OrthoFinder"]
'''

'''
# for BC-set metrics
methods_bc = ["BC-set", 
           "Neighborhood", 
           "Tree-Majority", 
           "Tree-Standard", 
           "Tree-Whitelist", 
           "PTP", 
           "OrthoFinder"]

rows_single=["Recall", "Precision", "F1-score", "Jaccard"]
'''

# for all reference
rows_all_ref = ["Neighborhood", 
        "Tree-Majority", 
        "Tree-Standard", 
        "Tree-Whitelist", 
        "PTP", 
        "OrthoFinder"]
cols_all_ref = ["Cactus", "BC-set", "B-set"] 


# for pure cactus metrics

'''
# for plants
methods_cactus = ["Cactus",
            "Neighborhood-1:1-gm", 
            "Original-orthologs",
            "PTP",
            "OrthoFinder"]
'''


# for drosophila
'''
methods_cactus = ["Cactus", 
                  "Neighborhood-global-pairs-gm",
            "Neighborhood-1:1-gm", 
            "OrthoFinder"]
'''

#rows_single=["Recall", "Precision", "F1-score", "Jaccard"]

'''
# plants metrics for pure cactus
metrics_cactus = np.array([[1.000, 0.403, 0.291, 0.096, 0.280], 
           [1.000, 0.413, 0.635, 0.351, 0.718],
           [1.000, 0.408, 0.399, 0.151, 0.403],
           [1.000, 0.256, 0.249, 0.081, 0.252]])
'''
           
# drosophila metrics for pure cactus
'''
metrics_cactus = np.array([[1.000, 0.112, 0.886, 0.786], 
           [1.000, 0.621, 0.907, 0.914],
           [1.000, 0.191, 0.896, 0.845],
           [1.000, 0.106, 0.812, 0.732]])
'''           

# ===============
# Mammals metrics
#================


# metrics as array (obtained from running evaluation scripts)
# Pure cactus
'''
metrics_cactus = np.array([[1.000, 0.857, 0.720, 0.723, 0.723, 0.624, 0.753], 
           [1.000, 0.855, 0.822, 0.825, 0.825, 0.861, 0.905],
           [1.000, 0.856, 0.768, 0.771, 0.771, 0.724, 0.822],
           [1.000, 0.748, 0.623, 0.627, 0.627, 0.567, 0.698]])
'''
'''
# BC-set
metrics_bc = np.array([[1.000, 0.759, 0.640, 0.643, 0.643, 0.555, 0.669], 
           [1.000, 0.865, 0.834, 0.837, 0.837, 0.874, 0.919],
           [1.000, 0.809, 0.724, 0.727, 0.727, 0.678, 0.774],
           [1.000, 0.679, 0.568, 0.571, 0.571, 0.513, 0.632]])
'''
# all reference
recall_all_ref = np.array([[0.857, 0.759, 0.710],
                            [0.720, 0.640, 0.644],
                            [0.723, 0.643, 0.646],
                            [0.723, 0.643, 0.646],
                            [0.624, 0.555, 0.595],
                            [0.753, 0.669, 0.675]])

precision_all_ref = np.array([[0.846, 0.856, 0.335],
                            [0.822, 0.834, 0.351],
                            [0.825, 0.837, 0.352],
                            [0.825, 0.837, 0.352],
                            [0.861, 0.874, 0.392],
                            [0.905, 0.919, 0.388]])

f1_all_ref = np.array([[0.852, 0.805, 0.455],
                            [0.768, 0.724, 0.454],
                            [0.771, 0.727, 0.456],
                            [0.771, 0.727, 0.456],
                            [0.724, 0.678, 0.472],
                            [0.822, 0.774, 0.493]])

jaccard_all_ref = np.array([[0.742, 0.673, 0.294],
                            [0.623, 0.568, 0.294],
                            [0.627, 0.571, 0.295],
                            [0.627, 0.571, 0.295],
                            [0.567, 0.513, 0.309],
                            [0.698, 0.632, 0.327]])


metrics = [
    ("Recall", recall_all_ref),
    ("Precision", precision_all_ref),
    ("F1-score", f1_all_ref),
    ("Jaccard", jaccard_all_ref)
]


'''
for metric in metrics_cactus:
    for value in metric:
        int(value)

for metric in metrics_bc:
    for value in metric:
        int(value)
'''        
for metric in recall_all_ref:
    for value in metric:
        int(value)
for metric in precision_all_ref:
    for value in metric:
        int(value)        
for metric in f1_all_ref:
    for value in metric:
        int(value) 
for metric in jaccard_all_ref:
    for value in metric:
        int(value)     
   


'''
# Plot all vs Cactus

fig, ax = plt.subplots(figsize=(10, 4), constrained_layout=True)

im = ax.imshow(
    metrics_cactus,
    cmap="viridis",
    vmin=0,
    vmax=1
)


# Show all ticks and label them with the respective list entries
ax.set_xticks(range(len(methods_cactus)), labels=methods_cactus,
              rotation=45, ha="right")
ax.set_yticks(range(len(rows_single)), labels=rows_single)

# Add value annotations
for i in range(len(rows_single)):
    for j in range(len(methods_cactus)):
        ax.text(j, i, f"{metrics_cactus[i, j]:.3f}",
                ha="center", va="center", color="black")

cbar = fig.colorbar(im, ax=ax, orientation="vertical", shrink=0.8)
cbar.set_label("Similarity Score", rotation=270, labelpad=15)

tick_values = np.linspace(0, 1, 6)  # 0.0 -> 1.0
cbar.set_ticks(tick_values)
cbar.set_ticklabels([f"{v:.1f}" for v in tick_values])

ax.set_title("Agreement between the methods of 1:1 ortholog evaluation\n"
             "vs reference (Cactus), Mammals")

# for plants

ax.set_title("Agreement between the methods of 1:1 ortholog evaluation\n"
             "vs reference (Cactus), Plants")

# for drosophila
ax.set_title("Agreement between the methods of 1:1 ortholog evaluation\n"
             "vs reference (Cactus), Drosophila")



# for mammals
plt.savefig("results/evaluation/metrics_vs_cactus_only.png",
    bbox_inches="tight",  # crop outer whitespace
    dpi=300)
# for plants

plt.savefig("results/evaluation/plants/metrics_vs_cactus_only.png",
    bbox_inches="tight",  # crop outer whitespace
    dpi=300)
    
# for drosophila

plt.savefig("results/evaluation/drosophila/metrics_vs_cactus_only.png",
    bbox_inches="tight",  # crop outer whitespace
    dpi=300)

plt.show()
plt.close()
'''

'''
# Plot all vs BC-set
fig, ax = plt.subplots(figsize=(8, 4))
im = ax.imshow(metrics_bc, cmap='viridis', vmin=0, vmax=1)

# Show all ticks and label them with the respective list entries
ax.set_xticks(range(len(methods_bc)), labels=methods_bc,
              rotation=45, ha="right")
ax.set_yticks(range(len(rows_single)), labels=rows_single)

# Add value annotations
for i in range(len(rows_single)):
    for j in range(len(methods_bc)):
        ax.text(j, i, f"{metrics_bc[i, j]:.3f}",
                ha="center", va="center", color="black")

cbar = fig.colorbar(im, ax=ax, orientation="vertical", shrink=0.8)
cbar.set_label("Similarity Score", rotation=270, labelpad=15)

tick_values = np.linspace(0, 1, 6)  # 0.0 -> 1.0
cbar.set_ticks(tick_values)
cbar.set_ticklabels([f"{v:.1f}" for v in tick_values])

ax.set_title("Agreement between the methods of 1:1 ortholog evaluation vs reference (BC-set)")
fig.tight_layout()
plt.savefig("results/evaluation/metrics_vs_bc_only.png")
plt.show()
plt.close()
'''

# Plot overview of three references, all metrics in one picture


fig, axs = plt.subplots(2, 2, figsize=(13, 13))
axs = axs.flatten()
last_im = None

for ax, (title, data) in zip(axs, metrics):

    im = ax.imshow(data, cmap='viridis', vmin=0, vmax=1)
    last_im = im

    ax.set_xticks(range(len(cols_all_ref)))
    ax.set_xticklabels(cols_all_ref, rotation=45, ha="right")

    ax.set_yticks(range(len(rows_all_ref)))
    ax.set_yticklabels(rows_all_ref)

    ax.set_title(title, fontsize=14)

    # Add values
    for i in range(len(rows_all_ref)):
        for j in range(len(cols_all_ref)):
            ax.text(j, i, f"{data[i, j]:.3f}",
                    ha="center", va="center", color="black")

fig.suptitle(
    "Agreement between methods of 1:1 ortholog evaluation vs reference (Cactus, BC-set and B-set)",
    fontsize=16,
    y=0.98   # move title slightly up
)
plt.tight_layout(rect=[0, 0, 0.92, 0.92])
plt.savefig("results/evaluation/mammalia/cactus_bc_and_b_heatmap.png")
plt.show()
plt.close()
