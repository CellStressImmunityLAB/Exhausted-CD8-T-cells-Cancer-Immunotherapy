Setting up the environment (from scratch)
==========================================
1) Download the correct Anaconda for your OS: https://www.anaconda.com
2) Download the following packages, using anaconda (all available through conda-forge): conda install -c conda-forge packagename
Package versions for the article are:
	multiprocessing-on-dill	3.5.0a4
	- numpy		1.18.4
	- pandas	0.25.3
	- python	3.7.0
	- scikit-learn	0.22.1
	- scipy		1.5.0
	- xgboost	0.9

The script was tested on Linux (Ubuntu 20.04.1 LTS Focal Fossa), macOS (10.14.6 Mojave) with updated versions of the packages listed.
We have not noticed any compatibility issues  as of 02/02/2021. Please report any problems running the script.

Installation time
===================
Typically within 10 minutes, depending on the speed of the Internet connection.
Less if these packages or equivalents are already installed.

Running the script and expected output 
=======================================
To run:
1) Go to the official link of the "The immune landscape of cancer": https://gdc.cancer.gov/about-data/publications/panimmune (REF: http://www.cell.com/immunity/retrieve/pii/S1074761318301213). You can also find this link in the "Data" folder.
2) In the data folder, place the two files necessary to run the script:
	- EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv (gene expression data)
	- TCGA.Kallisto.fullIDs.cibersort.relative_cancerImmuneLandscape.tsv (cell fraction data)
3) Change your current work directory to the folder containing the scripts. (e.g. cd /Users/MyUserName/Desktop/LOOCV_cspecific for macOSX)
4) Run the LOOCV.py script with Python in your favourite IDE/terminal: python LOOCV.py
5) Output should appear as an Excel table (CorrelationSignature_Preds_Real_Pearson_Garg.xlsx) in the current work directory and list the Pearson correlations per cancer type. If everything processed properly, you should see something like this:

	CancerType	PearsonCorr	Signature
0	BLCA	0.806333587	Garg
1	BRCA	0.723037027	Garg
2	GBM	0.465721993	Garg
3	LUAD	0.803623408	Garg
4	OV	0.715864752	Garg
5	SKCM	0.861922498	Garg


