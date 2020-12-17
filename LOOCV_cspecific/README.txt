Running the script
===================
To run:
1) Make a folder named "Data" in the directory containing the scripts.
2) In the data folder, place the two files from the cancer immune landscape paper:
	- EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv (gene expression data)
	- TCGA.Kallisto.fullIDs.cibersort.relative_cancerImmuneLandscape.tsv (cell fraction data)
	Files can be found here: https://gdc.cancer.gov/about-data/publications/panimmune
3) Make sure you have the requested packages installed:
	- Multiprocessing
	- Numpy
	- Pandas
	- Scikit-learn
	- Scipy
	- XGBoost
	These packages are often bundled with scientific Python installations, such as anaconda. 
	Packages not installed can be easily found through pip or in conda-forge.
	All other dependencies should come bundled with a default Python installation.
4) Change your current work directory to the folder containing the scripts.
5) Execute the LOOCV.py script with Python in your favourite IDE/terminal: python LOOCV.py
6) Output should appear as an Excel table in the current work directory and list the Pearson correlations per cancer type.
