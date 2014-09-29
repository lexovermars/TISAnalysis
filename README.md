TIS Analysis
===========


The identification of the correct translation initiation sites (TISs) constitutes an important aspect of sequence-based genome analysis. We have formulated a reference-free method to score the TIS annotation quality. The method is based on a comparison of the observed and expected distribution of all TISs in a particular genome given prior gene-calling. 

Requirements
============
- Python
- NumPy
- SciPy
- matplotlib
- R (only for the PCA-based correction of TIS existing annotations)

Usage
=====

TIS Annotation Quality
--------------

	python assess_TIS_annotation.py -i Example_data/NC_000913.ptt -f Example_data/NC_000913.fna -o "MyOutputName"

TIS Annotator
--------------
