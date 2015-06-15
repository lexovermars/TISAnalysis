TIS Analysis
===========


The identification of the correct translation initiation sites (TISs) constitutes an important aspect of sequence-based genome analysis. We have formulated a reference-free method to score the TIS annotation quality. The method is based on a comparison of the observed and expected distribution of all TISs in a particular genome given prior gene-calling. 

Requirements
============
- Python
- NumPy (python module)
- SciPy (python module)
- matplotlib (python module)

In addition, the PCa based correction existing TIS annotation requires:

- R 
- MASS (R package)

Usage
=====

TIS Annotation Quality
--------------

	python assess_TIS_annotation.py -i Example_Data/NC_000913.ptt -f Example_Data/NC_000913.fna -o "MyOutputName"
	
Output 

- MyOutputName_distribution.png: Distribution plot of genome-wide alternative TISs.
- MyOutputName_correlation.txt: Table with species name, GC-percentage, #ORFs and spearman correlation score


TIS Annotator
--------------
	 python run_pca_pipeline.py -i Example_data/NC_000913.ptt -f Example_data/NC_000913.fna -o "MyOutputName"


Output
- MyOutputName_5_rows_iteration_matrix.txt: Matrix with locus tag, TIS label (e.g. upstream, annotated or downstream) and PCA scores (PCA1,PCA2,PCA3; over 10 PCA iterations) for 5 top scoring TISs for each ORF. 
- MyOutputName_adjusted_annotation.txt Matrix with locus tag, TIS label (e.g. upstream, annotated or downstream) and PCA scores (PCA1,PCA2,PCA3; over 10 PCA iterations) for best scoring TIS. 
- MyOutputName_matrix_discarded.txt  Matrix with locus tag, TIS label (e.g. upstream, annotated or downstream) and PCA scores (PCA1,PCA2,PCA3; over 10 PCA iterations) for all TISs that did end up in the top 5 TIS table.
