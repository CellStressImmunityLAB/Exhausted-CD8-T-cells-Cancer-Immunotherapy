import numpy as np
import os
import pandas as pd
from sklearn.model_selection import LeaveOneOut
import multiprocessing as mp
from Support import runcv
from scipy.stats import pearsonr
import random

# Response variable to test signature on
response = "T.cells.CD8"

# ################# #
# Code starts below #
# ################# #
# Set seeds for reproducibility
random.seed(770)
np.random.seed(770)

# Literature analysis: 16 gene signature
signature = ["TNF", "IL2", "IFNG", "CD28", "B3GAT1", "TIGIT", "HAVCR2", "LAG3", "PDCD1",
             "EOMES", "CTLA4", "ENTPD1", "TOX", "CD244", "CD8A", "CD8B"]

# Immune landscape gene expression dataset
dgex = pd.read_csv(os.path.join("Data", "EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv"), sep="\t")
for somecol in list(dgex)[1:]:
    tempcol = []
    for x in dgex[somecol]:
        try:
            tempcol.append(np.float(x))
        except ValueError:
            tempcol.append(0)
    dgex[somecol] = tempcol
dgex.fillna(0, inplace=True)
siggex = dgex.loc[[True if x.split("|")[0] in signature else False for x in dgex["gene_id"]]]
siggex = siggex.transpose()
siggex.columns = siggex.iloc[0]
siggex = siggex.iloc[1:]
siggex.reset_index(drop=False, inplace=True)
keepcols = list(siggex)[1:]

# Immune landscape inferred fractions
cibersort = pd.read_csv(os.path.join("Data", "TCGA.Kallisto.fullIDs.cibersort.relative_cancerImmuneLandscape.tsv"),
                        sep="\t")
cibersort["SampleID"] = [x.replace(".", "-") for x in cibersort["SampleID"]]
cibersort["T.cells.CD8"] = pd.to_numeric(cibersort["T.cells.CD8"], errors="raise")
cibersort["T.cells.CD8"].fillna(0, inplace=True)

mergeme = pd.merge(left=siggex, right=cibersort, left_on="index", right_on="SampleID", how="inner")
mergeme.drop(columns=["index"], axis=1, inplace=True)

loocv = LeaveOneOut()

for somecol in keepcols:
    mergeme[somecol] = pd.to_numeric(mergeme[somecol], errors="raise")
mergeme["T.cells.CD8"] = pd.to_numeric(mergeme["T.cells.CD8"], errors="raise")
mergeme["T.cells.CD8"].fillna(0, inplace=True)
preds = []
reals = []

corrlist = []
ctypes = ["BRCA", "BLCA", "SKCM", "OV", "GBM", "LUAD"]

for ctype in sorted(ctypes):
    cframe = mergeme.loc[mergeme["CancerType"] == ctype]
    if len(cframe) > 50:
        myloocv = LeaveOneOut()
        mypool = mp.Pool(processes=60)
        myrun = [mypool.apply_async(runcv, args=(x, cframe, keepcols, response)) for x in
                 myloocv.split(cframe, groups=cframe["CancerType"])]
        output = [p.get() for p in myrun]
        mypool.close()
        lister = [[list(x[0])[0], x[1].values[0], x[2][0]] for x in output]
        CtypeRealPreds = pd.DataFrame(lister, columns=["CancerType", "Real", "Preds"])
        corrlist.append([ctype, pearsonr(CtypeRealPreds["Real"], CtypeRealPreds["Preds"])[0], "Garg"])
        print("Finished", "Garg", ctype)

corrframe = pd.DataFrame(corrlist, columns=["CancerType", "PearsonCorr", "Signature"])
corrframe.to_excel("CorrelationSignature_Preds_Real_Pearson_Garg.xlsx")
