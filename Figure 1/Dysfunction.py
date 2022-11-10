import pandas as pd
import numpy as np
import os
import matplotlib as mpl
mpl.use("tkagg")
import matplotlib.pyplot as plt
import seaborn as sns

fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(10, 3))
rootdir = os.path.join("Melanoma_correlationHeatmap_data_1658334863709")
i=0
myorder = ['naive.tsv', 'transitional.tsv', 'dysfunctional.tsv'] #os.listdir(rootdir)
remapper = {'naive': "Naive", 'transitional': "Transitional",'dysfunctional': "Antigen-specific\ndysfunction"}
for afile in myorder:
    fname = afile.split(".")[0]
    indata = pd.read_csv(os.path.join(rootdir, afile), sep=" ")
    sns.heatmap(indata, ax=ax[i], cmap="Reds", xticklabels=True, yticklabels=True, cbar_kws={"label":"Spearman correlation"})
    ax[i].set_title(remapper[fname])
    i += 1
plt.tight_layout()
plt.savefig("melanoma_comparison.pdf")
plt.close()

lister = []
indata = pd.read_csv(os.path.join("Melanoma_correlationHeatmap_data_1658334863709",
                                  "Z-score_Mean_No. cells_with ZEROES_data_Melanoma.tsv"), sep="\t")
indata = indata.groupby(by="x").mean()

biodic = {}
backmap = {}
biomartdata = pd.read_csv("mart_export_11102021.txt", sep="\t")
for x, y in zip(biomartdata["Gene name"], biomartdata["Gene stable ID"]):
    biodic[str(x)] = str(y)
    backmap[str(y)] = str(x)

renamer = {}
for acol in list(indata):
    renamer[acol] = backmap[acol.replace("gene_", "")]
indata.rename(columns=renamer, inplace=True)

barframe = indata.copy()
for acol in list(barframe):
    mysum = np.sum(barframe[acol])
    barframe[acol] = [x / mysum for x in list(barframe[acol])]
barframe = barframe.transpose()
barframe.reset_index(drop=False, inplace=True)
barframe.sort_values(by="dysfunctional", inplace=True)
barframe = barframe[["index", "naive", "transitional", "dysfunctional"]].copy()
barframe.rename(columns={"naive": "Naive", "transitional": "Transitional",
                         "dysfunctional": "Antigen-specific\ndysfunction"}, inplace=True)
sigmapper = {"Naive": sns.color_palette()[0], "Transitional": sns.color_palette()[1],
             "Antigen-specific\ndysfunction": sns.color_palette()[2]}
fig, ax = plt.subplots()
i = 0
prev = [0]
for avar in sorted(set(list(barframe)[1:])):
    if i == 0:
        ax.bar(list(barframe["index"]), list(barframe[avar]), color=sigmapper[avar], label=avar)
        prev = [0] * len(barframe)
    else:
        ax.bar(list(barframe["index"]), list(barframe[avar]), bottom=prev, color=sigmapper[avar], label=avar)
    prev = prev + barframe[avar]
    i = i + 1
ax.set_ylabel('Proportional median expression')
ax.legend(bbox_to_anchor=(1.03, 1))
plt.xticks(rotation=90)
plt.tight_layout()
sns.despine()
plt.savefig("Stacked_barplot_genes_dysfunction.pdf")
plt.close()
