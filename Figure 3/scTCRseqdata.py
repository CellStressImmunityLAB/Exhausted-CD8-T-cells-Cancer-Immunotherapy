import itertools
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("TkAgg")
import scanpy as sc
import scirpy as ir
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact
from sklearn.preprocessing import robust_scale
from adjustText import adjust_text
from collections import Counter
import scipy
from statannot import add_stat_annotation
from scipy.stats import mannwhitneyu
vdjdb = ir.datasets.vdjdb()

# Glioma
samples = os.listdir(os.path.join("Data", "TCRseq", "GBM", "GSE163108_RAW"))
i = 0
for asample in samples:
    prefix = asample.split("_")[0]+"_GEX_"
    exprdata = sc.read_10x_h5(os.path.join("Data", "TCRseq", "GBM", "GSE163108_RAW", asample))
    exprdata.var_names_make_unique()
    exprdata.obs.index = [prefix+x[:-2] for x in exprdata.obs.index]
    if i == 0:
        glioma = exprdata.copy()
    else:
        glioma = sc.concat([glioma, exprdata])
    i = i + 1

vdjdata = pd.read_csv(os.path.join("Data", "TCRseq", "GBM", "GSE163108_metadata_10x.csv"), sep=",")
vdjdata = vdjdata.loc[vdjdata["annotate_Tcelltype"] == "CD8+"].copy()
vdjdata.dropna(subset=["cdr3s_aa"], inplace=True)
vdjdata["numcd3"] = [len(x.split(";")) for x in vdjdata["cdr3s_aa"]]
vdjdata = vdjdata.loc[vdjdata["numcd3"] < 3].copy()
vdjdata = vdjdata.loc[vdjdata["annotate_clonotype"] == "CD8+"].copy()
vdjdata = vdjdata.loc[~vdjdata["sampleid"].isin(["E37_GEX", "E99_GEX", "E100_GEX"])].copy()
glioma = glioma[sorted(set(vdjdata["Unnamed: 0"])), :].copy()
vdjdata.set_index("Unnamed: 0", inplace=True)

tcr_cells = []
for idx, row in vdjdata.iterrows():
    cell = ir.io.AirrCell(cell_id=idx)
    alpha_chain = ir.io.AirrCell.empty_chain_dict()
    beta_chain = ir.io.AirrCell.empty_chain_dict()

    elements = row["cdr3s_aa"].split(";")
    for element in elements:
        if element.startswith("TRA:"):
            print(element)
            d = {"locus": "TRA",
                 "junction_aa": element.replace("TRA:", ""),
                 "productive": "True"}
            alpha_chain.update(d)
            cell.add_chain(alpha_chain)
        else:
            if element.startswith("TRB:"):
                e = {"locus": "TRB",
                     "junction_aa": element.replace("TRB:", ""),
                     "productive": "True"}
                beta_chain.update(e)
                cell.add_chain(beta_chain)
            else:
                pass
    tcr_cells.append(cell)
    try:
        del d
    except NameError:
        pass
    try:
        del e
    except NameError:
        pass

vdjdata = ir.io.from_airr_cells(tcr_cells)
ir.pp.merge_with_ir(glioma, vdjdata)
glioma.obs["sample"] = [x.split("_")[0] for x in list(glioma.obs.index)]
glioma.obs["Cancer"] = "GBM"

# BRCAdata
samples = os.listdir(os.path.join("Data", "TCRseq", "BRCA", "GSE114724_RAW"))
i = 0
for asample in samples:
    exprdata = sc.read_10x_mtx(os.path.join("Data", "TCRseq", "BRCA", "GSE114724_RAW", asample))
    vdjdata = ir.io.read_10x_vdj(os.path.join("Data", "TCRseq", "BRCA", "GSE114724_RAW", asample,
                                              "filtered_contig_annotations.csv.gz"))
    ir.pp.merge_with_ir(exprdata, vdjdata)
    exprdata.obs["sample"] = asample
    if i == 0:
        brcadata = exprdata.copy()
    else:
        brcadata = sc.concat([brcadata, exprdata])
    i = i+1
brcadata.obs["Cancer"] = "BRCA"

# LUADdata
luaddata = sc.read_csv(os.path.join("Data", "TCRseq", "LUAD", "tcellcounts.csv"))
luaddata = luaddata.transpose()
luaddata.obs["sample"] = ["-".join(x.split(".")[:2]) for x in list(luaddata.obs.index)]

vdjdata = pd.read_csv(os.path.join("Data", "TCRseq", "LUAD", "GSE179994_all.scTCR.tsv"), sep="\t")
vdjdata.set_index("CellName", inplace=True)
vdjdata["sample"] = ["-".join(x.split(".")[:2]) for x in list(vdjdata.index)]

# Keep only untreated
keepsamps = [x for x in set(luaddata.obs["sample"]) if x.find("-ut") != -1]
keepsamps = [x for x in set(vdjdata["sample"]) if x in keepsamps]
vdjdata = vdjdata.loc[vdjdata["sample"].isin(keepsamps)].copy()
keepindexes = set(luaddata.obs.index).intersection(set(vdjdata.index))
vdjdata = vdjdata.loc[list(keepindexes), :].copy()
luaddata = luaddata[list(keepindexes), :].copy()

# Create a list of AnnData objects (one for each sample)
tcr_cells = []
for idx, row in vdjdata.iterrows():
    cell = ir.io.AirrCell(cell_id=idx)
    alpha_chain = ir.io.AirrCell.empty_chain_dict()
    beta_chain = ir.io.AirrCell.empty_chain_dict()

    d = {
        "locus": "TRA",
        "junction_aa": row["CDR3(Alpha1)"],
        "junction": row["CDR3_nt(Alpha1)"],
        "consensus_count": float(row["nRead(Alpha1)"]),
        "v_call": row["V_gene(Alpha1)"],
        # "d_call": "None",
        "j_call": row["J_gene(Alpha1)"],
        "productive": "True",
    }
    alpha_chain.update(d)
    cell.add_chain(alpha_chain)

    d = {
        "locus": "TRB",
        "junction_aa": row["CDR3(Beta1)"],
        "junction": row["CDR3_nt(Beta1)"],
        "consensus_count": float(row["nRead(Beta1)"]),
        "v_call": row["V_gene(Beta1)"],
        "j_call": row["J_gene(Beta1)"],
        "productive": "True",
    }
    beta_chain.update(d)
    cell.add_chain(beta_chain)
    tcr_cells.append(cell)

vdjdata = ir.io.from_airr_cells(tcr_cells)
ir.pp.merge_with_ir(luaddata, vdjdata)
luaddata.obs["Cancer"] = "LUAD"

# SKCMdata
samples = os.listdir(os.path.join("Data", "TCRseq", "SKCM", "Mergedata"))
i = 0
for asample in samples:
    exprdata = sc.read_10x_mtx(os.path.join("Data", "TCRseq", "SKCM", "Mergedata", asample))
    vdjdata = ir.io.read_10x_vdj(os.path.join("Data", "TCRseq", "SKCM", "Mergedata", asample,
                                              "filtered_contig_annotations.csv.gz"))
    ir.pp.merge_with_ir(exprdata, vdjdata)
    exprdata.obs["sample"] = asample

    if i == 0:
        skcmdata = exprdata.copy()
    else:
        skcmdata = sc.concat([skcmdata, exprdata])
    i = i + 1

for afile in os.listdir(os.path.join("Data", "TCRseq", "SKCM", "GSE159251_RAW")):
    if not afile.startswith("."):
        exprdata = sc.read_h5ad(os.path.join("Data", "TCRseq", "SKCM", "GSE159251_RAW", afile, afile+".h5ad"))
        vdjdata = pd.read_csv(os.path.join("Data", "TCRseq", "SKCM", "GSE159251_RAW", afile, afile+".txt"), sep=" ")
        exprdata.obs["sample"] = afile
        # Create a list of AnnData objects (one for each sample)
        tcr_cells = []
        for idx, row in vdjdata.iterrows():
            cell = ir.io.AirrCell(cell_id=idx)
            alpha_chain = ir.io.AirrCell.empty_chain_dict()
            beta_chain = ir.io.AirrCell.empty_chain_dict()

            if row["Matching"] == "notcr":
                cell.add_chain(alpha_chain)
                cell.add_chain(beta_chain)
            else:
                if row["Matching"] == "alpha_chain_only":
                    d = {
                        "locus": "TRA",
                        "junction_aa": row["TCR"].replace("|", ""),
                        "consensus_count": float(row["nCounts"]),
                        "productive": "True"
                    }
                    alpha_chain.update(d)
                    cell.add_chain(alpha_chain)
                else:
                    if row["Matching"] == "beta_chain_only":
                        d = {
                            "locus": "TRB",
                            "junction_aa": row["TCR"].replace("|", ""),
                            "consensus_count": float(row["nCounts"]),
                            "productive": "True"
                        }
                        beta_chain.update(d)
                        cell.add_chain(beta_chain)
                    else:
                        d = {
                            "locus": "TRA",
                            "junction_aa": row["TCR"].split("|")[0],
                            "consensus_count": float(row["nCounts"]),
                            "productive": "True",
                        }
                        alpha_chain.update(d)
                        cell.add_chain(alpha_chain)

                        d = {
                            "locus": "TRB",
                            "junction_aa": row["TCR"].split("|")[1],
                            "consensus_count": float(row["nCounts"]),
                            "productive": "True",
                        }
                        beta_chain.update(d)
                        cell.add_chain(beta_chain)
            tcr_cells.append(cell)

        vdjdata = ir.io.from_airr_cells(tcr_cells)
        ir.pp.merge_with_ir(exprdata, vdjdata)
        skcmdata = sc.concat([skcmdata, exprdata])
skcmdata.obs["Cancer"] = "SKCM"
processdata = sc.concat([glioma, brcadata, luaddata, skcmdata])

processdata.obs_names_make_unique()
sc.pp.calculate_qc_metrics(processdata, inplace=True)
sc.pp.filter_cells(processdata, min_genes=250, inplace=True)
sc.pp.filter_cells(processdata, max_genes=2500, inplace=True)
sc.pp.filter_genes(processdata, min_cells=1000, inplace=True)
sc.pp.normalize_per_cell(processdata, counts_per_cell_after=10000)
sc.pp.log1p(processdata)

keepindexes = processdata.obs.loc[processdata.obs["has_ir"] == "True"].index
processdata = processdata[keepindexes, :].copy()

tempdf = processdata[:, ["CD8A", "CD8B"]].to_df()
tempdf["sum"] = tempdf.sum(axis=1)
tempdf = tempdf.loc[tempdf["sum"] > 0].copy()
keepindexes = tempdf.index
processdata = processdata[keepindexes, :].copy()

sc.pp.combat(processdata, key="sample")

# csCD metagene
sigtable = pd.read_excel(os.path.join("Data", "TableX_signature_TCRseq.xlsx"), engine="openpyxl", skiprows=2)
sigdic = {}
for acol in list(sigtable):
    lister = []
    for agene in set(sigtable[acol]):
        if not pd.isna(agene):
            if agene in processdata.var_names:
                lister.append(agene)
            else:
                pass
    if len(lister) > 1:
        sigdic[acol] = lister

for asig in sigdic.keys():
    getter = processdata[:, sigdic[asig]].to_df()
    processdata.obs[asig] = getter.sum(axis=1)/len(list(getter))

sc.pp.highly_variable_genes(processdata, flavor="cell_ranger", n_top_genes=2000)
sc.tl.pca(processdata)
sc.pp.neighbors(processdata)
sc.tl.leiden(processdata)
sc.tl.umap(processdata, n_components=2)
sc.pl.umap(processdata, color=["leiden", "Cancer"])
sns.despine()
plt.autoscale()
plt.savefig("UMAP_cancer.png")
plt.close()

ir.tl.chain_qc(processdata)
keepindexes = processdata.obs.loc[processdata.obs["multi_chain"] == "False"].index
processdata = processdata[keepindexes, :].copy()

ax = ir.pl.group_abundance(processdata, groupby="chain_pairing", target_col="Cancer")
sns.despine()
plt.tight_layout()
plt.tight_layout()
plt.savefig("chain_pairing_count.pdf")
plt.close()

# TRA + TRB
vdic = {}
i = 0
lister = []
for x, y in zip(processdata.obs["IR_VJ_1_junction_aa"], processdata.obs["IR_VDJ_1_junction_aa"]):
    if not pd.isna(x):
        try:
            lister.append(vdic[x])
        except KeyError:
            vdic[x] = str(i)
            if not pd.isna(y):
                vdic[y] = str(i)
            lister.append(str(i))
            i = i+1
    else:
        try:
            lister.append(vdic[y])
        except KeyError:
            vdic[y] = str(i)
            lister.append(str(i))
            i = i+1

processdata.obs["clone_id"] = lister

# Define thresholding
ir.tl.clonotype_modularity(processdata, target_col='clone_id', connectivity_key='connectivities',
                           permutation_test='approx', n_permutations=None, key_added='clonotype_modularity',
                           inplace=True, fdr_correction=True, random_state=0)

lister = []
expmod = []
for athres in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
    ir.tl.clonal_expansion(processdata, expanded_in="Cancer", target_col="clone_id", clip_at=athres, inplace=True,
                           key_added="Expansion")
    processdata.obs["Expansion"] = ["Expander" if x == ">= "+str(athres) else "Non-expander" for x in
                                    list(processdata.obs["Expansion"])]
    expmod.append(np.mean(processdata.obs["clonotype_modularity"].loc[processdata.obs["Expansion"] == "Expander"]))
    for acancer in set(processdata.obs["Cancer"]):
        mycancer = processdata.obs.loc[processdata.obs["Cancer"] == acancer].copy()
        mycounts = Counter(mycancer["Expansion"])
        totals = sum(mycounts.values())
        lister.append([acancer, athres, "Non-expander", mycounts["Non-expander"], mycounts["Non-expander"]/totals])
        lister.append([acancer, athres, "Expander", mycounts["Expander"], mycounts["Expander"]/totals])
thresframe = pd.DataFrame(lister, columns=["Cancer", "Threshold", "Type", "Count", "Proportion"])

pframe = thresframe.loc[thresframe["Type"] == "Expander"].copy()

# Plot expanders vs non-expanders
gets = pframe.pivot_table(index="Threshold", columns="Cancer", values="Proportion")
gets["All"] = (gets["BRCA"]+gets["LUAD"]+gets["SKCM"])/3
gets = gets[["GBM", "All"]].copy()

yvals = [int(x) for x in gets.index]
gbmvals = gets["GBM"]

gets = pframe.pivot_table(index="Threshold", columns="Cancer", values="Proportion")
yvals = gets.index
plt.scatter(gets["GBM"].values, yvals, color=sns.color_palette()[0], alpha=1, label='GBM')

lister = []
for athres in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
    ir.tl.clonal_expansion(processdata, expanded_in="Cancer", target_col="clone_id", clip_at=athres,
                           inplace=True, key_added="Expansion")
    processdata.obs["Expansion"] = ["Expander" if x == ">= "+str(athres) else "Non-expander" for x
                                    in list(processdata.obs["Expansion"])]

    gbmmer = processdata.obs.loc[processdata.obs["Cancer"] == "GBM"]
    other = processdata.obs.loc[processdata.obs["Cancer"] != "GBM"]

    gbmmer = Counter(gbmmer["Expansion"])
    other = Counter(other["Expansion"])

    lister.append([str(athres), other["Expander"]/(other["Expander"]+other["Non-expander"]),
                   fisher_exact([[other["Expander"], other["Non-expander"]], [gbmmer["Expander"],
                                                                              gbmmer["Non-expander"]]])[1]])

lframe = pd.DataFrame(lister, columns=["Threshold", "Contrast", "Pval"])
lframe["Pval"] = ["*p="+str(np.format_float_scientific(pval, 2)) if pval < 0.05 else
                  "p="+str(np.format_float_scientific(pval, 2)) for pval in lframe["Pval"]]
lframe.to_csv("pvals_thres.tsv", sep="\t")

i = 0
for xval, yval, pval in zip(lframe["Contrast"], lframe["Threshold"], lframe["Pval"]):
    if int(yval) > 1:
        try:
            plt.text(xval+0.01, float(yval)+0.2, s=pval)
        except IndexError:
            plt.text(xval+0.01, float(yval)+0.2, s=pval)
        i = i + 1

i = 0
for acol in ["BRCA", "LUAD", "SKCM"]:
    if not acol == "GBM":
        plt.scatter(gets[acol].values, yvals, color=sns.color_palette()[i], alpha=0.4, label=acol)
    i = i + 1
plt.hlines(y=yvals, xmin=gets["GBM"].values,
           xmax=[max(x[0], x[1]) for x in zip(gets["LUAD"].values, gets["SKCM"].values)], color='grey', alpha=0.4)

plt.yticks(np.arange(0, 11))
plt.legend()
plt.ylim(1.5, 11)
plt.xlim(0.38, 0.75)
plt.tight_layout()
sns.despine()
plt.savefig("thresholdselection_TCR.pdf")
plt.close()

# 5 as the chosen threshold
athres = 5
ir.tl.clonal_expansion(processdata, expanded_in="Cancer", target_col="clone_id", clip_at=athres, inplace=True,
                       key_added="Expansion")
processdata.obs["Expansion"] = ["Expander" if x == ">= " + str(athres) else "Non-expander" for x in
                                list(processdata.obs["Expansion"])]
fig, ax = plt.subplots()
xvals = []
yvals = [[], []]
i = 0
for acancer in ["GBM", "BRCA", "LUAD", "SKCM"]:
    subframe = processdata.obs.loc[processdata.obs["Cancer"] == acancer]
    mycounts = Counter(subframe["Expansion"])
    print(mycounts)
    xvals.append(acancer)
    yvals[0].append(mycounts["Expander"]/(mycounts["Expander"]+mycounts["Non-expander"]))
    yvals[1].append(mycounts["Non-expander"]/(mycounts["Expander"]+mycounts["Non-expander"]))
    i = i + 1
plt.stackplot(xvals, yvals, labels=["Expander", "Non-expander"])
plt.legend(loc='upper left')

i = 0
for acancer in ["GBM", "BRCA", "LUAD", "SKCM"]:
    plt.text(acancer, yvals[0][i], str(np.round(yvals[0][i], 3)))
    i = i + 1
plt.tight_layout()
sns.despine()
plt.savefig("ExpanderProps.pdf")
plt.close()

# Expansion umap
processdata.obs["UMAP1"] = [x[0] for x in processdata.obsm["X_umap"]]
processdata.obs["UMAP2"] = [x[1] for x in processdata.obsm["X_umap"]]
fig, ax = plt.subplots(figsize=(4, 3))
sns.kdeplot(x="UMAP1", y="UMAP2", data=processdata.obs, hue="Expansion", ax=ax, fill=True,
            thresh=0.60, levels=[0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1],
            hue_order=["Expander", "Non-expander"],
            palette={"Expander": sns.color_palette()[9], "Non-expander": sns.color_palette()[7]})
plt.tight_layout()
sns.despine()
plt.savefig("kdeplot_expansion.pdf")
plt.close()

fig, ax = plt.subplots(figsize=(4, 3))
sc.pl.umap(processdata, color=["Expansion"], palette={"Expander": sns.color_palette()[9],
                                                      "Non-expander": sns.color_palette()[7]}, ax=ax)
sns.despine()
plt.tight_layout()
plt.savefig("umap_expansion.png")
plt.close()

fig, ax = plt.subplots(figsize=(6, 3))
sc.pl.umap(processdata, color=["Cancer"], palette={"GBM": sns.color_palette()[0],
                                                   "SKCM": sns.color_palette()[1], "LUAD": sns.color_palette()[2],
                                                   "BRCA": sns.color_palette()[3]}, ax=ax)
sns.despine()
plt.tight_layout()
plt.savefig("umap_cancer.png")
plt.close()

fig, ax = plt.subplots(figsize=(6, 3))
sc.pl.umap(processdata, color=["sample"], ax=ax)
sns.despine()
plt.tight_layout()
plt.savefig("umap_sample.png")
plt.close()

for asig in ['IFNg score', 'Expanded immune gene signature', 'T cell–inflamed GEP', 'IMPRES signature', 'CYT score',
             'T exhaust signature', 'T accum signature', 'CTL signature', 'TIDE signature', 'csCD8 signature',
             'Differential_exhaustion', 'NeoTCR8', 'Bystander_activated', 'TCR_activated', 'Severe_dysfunction']:

    fig, ax = plt.subplots(figsize=(3, 4))
    sns.violinplot(x="Cancer", y=asig, hue="Expansion", data=processdata.obs, order=["SKCM", "LUAD", "BRCA", "GBM"],
                   hue_order=["Expander", "Non-expander"], clip=0, split=True,
                   palette={"Expander": sns.color_palette()[1], "Non-expander": sns.color_palette()[0]},
                   inner="quartile", ax=ax)
    pvals = []
    boxpairs = []
    for acancer in ["SKCM", "LUAD", "BRCA", "GBM"]:
        subframe = processdata.obs.loc[processdata.obs["Cancer"] == acancer].copy()
        g1 = subframe.loc[subframe["Expansion"] == "Expander"]
        g2 = subframe.loc[subframe["Expansion"] == "Non-expander"]
        log2fc = np.log2(np.mean(g1[asig])) - np.log2(np.mean(g2[asig]))
        if log2fc > 0.1:
            boxpairs.append(((acancer, "Expander"), (acancer, "Non-expander")))
            pvals.append(mannwhitneyu(g1[asig], g2[asig])[1])
        print(acancer, len(subframe), log2fc)

    g1 = processdata.obs.loc[(processdata.obs["Expansion"] == "Expander") & (processdata.obs["Cancer"] == "GBM")]
    g2 = processdata.obs.loc[(processdata.obs["Expansion"] == "Expander") & (processdata.obs["Cancer"] == "LUAD")]
    pvals.append(mannwhitneyu(g1[asig], g2[asig])[1])
    boxpairs.append((("LUAD", "Expander"), ("GBM", "Expander")))
    g1 = processdata.obs.loc[(processdata.obs["Expansion"] == "Expander") & (processdata.obs["Cancer"] == "GBM")]
    g2 = processdata.obs.loc[(processdata.obs["Expansion"] == "Expander") & (processdata.obs["Cancer"] == "SKCM")]
    pvals.append(mannwhitneyu(g1[asig], g2[asig])[1])
    boxpairs.append((("SKCM", "Expander"), ("GBM", "Expander")))

    add_stat_annotation(ax=ax, x="Cancer", y=asig, hue="Expansion", data=processdata.obs,
                        order=["SKCM", "LUAD", "BRCA", "GBM"], hue_order=["Expander", "Non-expander"],
                        perform_stat_test=False, text_format="full", pvalues=pvals,
                        box_pairs=boxpairs)
    plt.xlabel("")
    plt.tight_layout()
    sns.despine()
    plt.savefig("violin_expanders_"+asig.replace("/", "-")+"cutofflog2fc0point1.pdf")
    plt.close()

clonomap = {}
for clontype in set(processdata.obs["clone_id"]):
    subf = processdata.obs.loc[processdata.obs["clone_id"] == clontype].copy()
    thecounts = Counter(subf["Cancer"])
    thecombos = [x for x in itertools.combinations(thecounts.keys(), 2)]

    if len(thecounts.keys()) == 1:
        try:
            clonomap[(list(thecounts.keys())[0], list(thecounts.keys())[0])] = clonomap[(list(thecounts.keys())[0],
                                                                                         list(thecounts.keys())[0])]+1
        except KeyError:
            clonomap[(list(thecounts.keys())[0], list(thecounts.keys())[0])] = 1
    else:
        for acombo in thecombos:
            try:
                clonomap[acombo] = clonomap[acombo]+1
            except KeyError:
                clonomap[acombo] = 1

# Chord diagram
lister = []
for akey in clonomap.keys():
    lister.append([akey[0], akey[1], clonomap[akey]])
thedf = pd.DataFrame(lister, columns=["Source", "Target", "Value"])
thedf.to_csv("cytoscape_interactions.txt", sep="\t")

sc.tl.rank_genes_groups(processdata, groupby="Expansion")

selected = pd.read_excel(os.path.join("Data", "genelist_N_v_NE_ADG.xlsx"), engine="openpyxl")
selected = set(selected["Gene"].loc[selected["Selected"] == 1])
group0 = []
mylabels = []
colorder = processdata.uns["rank_genes_groups"]["logfoldchanges"].dtype.names
for x, y, z, a in zip(processdata.uns["rank_genes_groups"]["names"],
                      processdata.uns["rank_genes_groups"]["logfoldchanges"],
                      processdata.uns["rank_genes_groups"]["pvals_adj"],
                      processdata.uns["rank_genes_groups"]["scores"]):
    mycolor = ""
    if y[0] > 0:
        mycolor = "crimson"
    else:
        if y[0] < 0:
            mycolor = "teal"
        else:
            mycolor = "grey"
    if np.abs(y[0]) > 0.25:
        if x[0] in selected:
            group0.append([x[0], y[0], z[0], a[0], x[0], mycolor])
        else:
            group0.append([x[0], y[0], z[0], a[0], "", mycolor])
    else:
        group0.append([x[0], y[0], z[0], a[0], "", mycolor])

group0 = pd.DataFrame(group0, columns=["Gene", "LogFC", "Padj", "Score", "Label", "Color"])
group0["Group"] = colorder[0]
group0["Padj"] = -group0["Padj"].apply(np.log10)
themax = max([x for x in list(group0["Padj"]) if np.isfinite(x)])
group0["Padj"] = [x if np.isfinite(x) else themax for x in list(group0["Padj"])]
group0.rename(columns={"Padj": "-log10(adjusted P-value)"}, inplace=True)

fig, ax = plt.subplots(figsize=(5, 3))
xvals = group0["LogFC"]
yvals = group0["-log10(adjusted P-value)"]
ax.scatter(x=xvals, y=yvals, c=group0["Color"], s=8)
texts = []
for x, y, z in zip(xvals, yvals, group0["Label"]):
    if not z == "":
        texts.append(plt.text(x, y, z, size=8))
plt.axvline(0, linewidth=0.2)
plt.axhline(-np.log10(0.05), linestyle="dashed", color="black", alpha=0.6)
plt.xlabel("log2(fold change)")
plt.ylabel("-log10(adjusted P-value)")
adjust_text(texts, arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
plt.tight_layout()
sns.despine()
plt.savefig("volcanoplot_tcr.pdf")
plt.close()

gets = group0.copy()
gets = gets.loc[gets["Gene"].isin(processdata.var["highly_variable"].loc[
                                      processdata.var["highly_variable"]].index)].copy()
gets.reset_index(drop=True, inplace=True)

panimmune = pd.read_excel(os.path.join("Data", "PanImmune_GeneSet_Definitions.xlsx"),
                          engine="openpyxl", sheet_name="Genes")
thegenes = set(panimmune["Gene"])

gets = gets.loc[gets["Gene"].isin(thegenes)].copy()
gets.reset_index(drop=True, inplace=True)

expandersig = gets["Gene"].iloc[:30].copy()
nonexpandersig = gets["Gene"].iloc[-30:].copy()

fx = open("expanders_gmt.gmt", "w")
fx.write("Expanders\tExpanders\t" + "\t".join(sorted(set(expandersig)))+"\n")
fx.write("Non-expanders\tNon-expanders\t" + "\t".join(sorted(set(nonexpandersig))))
fx.close()

expanders = processdata[:, [x for x in set(expandersig)]].to_df()
expanders = pd.DataFrame(expanders.mean(axis=1), columns=["E"])
processdata.obs["E"] = list(expanders["E"])

nonexpanders = processdata[:, [x for x in set(nonexpandersig)]].to_df()
nonexpanders = pd.DataFrame(nonexpanders.mean(axis=1), columns=["NE"])
processdata.obs["NE"] = list(nonexpanders["NE"])

lister = []
for asig in ["SKCM", "LUAD", "BRCA", "GBM"]:
    subf = processdata.obs.loc[(processdata.obs["Cancer"] == asig) & (processdata.obs["Expansion"] == "Expander")]
    evals = [x for x in subf["E"] if x > 0]
    myrange = scipy.stats.iqr(evals)*1.5
    top = np.quantile(evals, 0.75)
    evals = [x for x in evals if x < (top+myrange)]
    nevals = [x for x in subf["NE"] if x > 0]
    myrange = scipy.stats.iqr(nevals)*1.5
    top = np.quantile(nevals, 0.75)
    nevals = [x for x in nevals if x < (top+myrange)]
    lister.append([asig, np.median(evals), np.median(nevals)])

lframe = pd.DataFrame(lister, columns=["Cancer", "E", "NE"])
lframe.set_index("Cancer", inplace=True)

lframe = pd.DataFrame(robust_scale(lframe, axis=0), columns=list(lframe), index=lframe.index)
sns.heatmap(lframe, cmap="RdYlBu_r", cbar_kws={"label": "IQR-normalised signature expression"})
plt.savefig("robust_scaled.pdf")
plt.close()
