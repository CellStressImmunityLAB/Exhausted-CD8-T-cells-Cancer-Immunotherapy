import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns

rootdir = os.path.join("HIV_correlationHeatmap_data_1658337517096")
i = 0
myorder = ['0.3 years.tsv', '3.5 years.tsv', '7.6 years.tsv', '11 years.tsv']
remapper={'0.3 years':'CD8+T (0.3 years)', '3.5 years':'CD8+T (3.5 years)',
          '7.6 years':'CD8+T (7.6 years)', '11 years':'CD8+T (11 years)'}

fig, ax = plt.subplots(nrows=1, ncols=4, figsize=(14, 3))
for afile in myorder:
    fname=afile.replace(".tsv","")
    indata=pd.read_csv(os.path.join(rootdir, afile), sep=" ")
    sns.heatmap(indata, ax=ax[i], cmap="Reds", xticklabels=True, yticklabels=True, cbar_kws={"label":"Spearman correlation"})
    ax[i].set_title(remapper[fname])
    i+=1
plt.tight_layout()
plt.savefig("chronic HIV.pdf")
plt.close()

lister=[]
indata=pd.read_csv(os.path.join(rootdir, "Z-score_Mean_No. cells_with ZEROES_data.tsv"), sep="\t")
indata=indata.groupby(by="x").mean()

biodic = {}
backmap = {}
biomartdata = pd.read_csv("mart_export_11102021.txt", sep="\t")
for x, y in zip(biomartdata["Gene name"], biomartdata["Gene stable ID"]):
    biodic[str(x)] = str(y)
    backmap[str(y)] = str(x)

renamer={}
for acol in list(indata):
    renamer[acol]=backmap[acol.replace("gene_","")]
indata.rename(columns=renamer, inplace=True)

barframe=indata.copy()
for acol in list(barframe):
    mysum = np.sum(barframe[acol])
    barframe[acol] = [x / mysum for x in list(barframe[acol])]

barframe = barframe.transpose()
barframe.reset_index(drop=False, inplace=True)
barframe.rename(columns=remapper, inplace=True)
barframe.sort_values(by="CD8+T (3.5 years)", inplace=True)

barframe=barframe[["index", "CD8+T (0.3 years)", "CD8+T (3.5 years)", "CD8+T (7.6 years)", "CD8+T (11 years)"]].copy()
sigmapper = {"CD8+T (0.3 years)": sns.color_palette()[0], "CD8+T (3.5 years)": sns.color_palette()[1],
             "CD8+T (7.6 years)": sns.color_palette()[2], "CD8+T (11 years)":sns.color_palette()[3]}

fig, ax = plt.subplots()
i = 0
prev = [0]
for avar in list(barframe)[1:]:
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
plt.savefig("Stacked_barplot_genes_HIV.pdf")
plt.close()
