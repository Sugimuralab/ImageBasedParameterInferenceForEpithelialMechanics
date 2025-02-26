Parameter estimation and model selection for epithelial mechanics
===

## Description

Image-basde parameter inference has been reported in Ogita et al. 2022 [1]. It estimates parameters of the Cell Vertex Model (CVM) and selects the most predictive model among tested from image data of epithelial tissue.
See Ogita et al. 2022 [1] for details.

The script for Bayesian force inference can be found [here](https://github.com/IshiharaLab/BayesianForceInference).

Using a input file that contains the information about the position and connectivity of cell vertices from an image of epithelial tissue, 
the scripts can be used to perform image-based parameter inference with "ParameterEstimation.py". 


## Requirement

* statsmodels 0.12.2


## Usage

1. Prepare an input file in the same format as the attached sample (Sample/sample/).
2. Change the variable "filename" in ParameterEstimation.py to the input file in step 1.
3. Run ParameterEstimation.py on IDE or IPython.


## Reference

1. Goshi Ogita, Takefumi Kondo, Keisuke Ikawa, Tadashi Uemura, Shuji Ishihara & Kaoru Sugimura (2022)<br>
"Image based parameter inference for epithelial mechanics"<br>
PLoS Comput Biol 18(6): e1010209. [https://doi.org/10.1371/journal.pcbi.1010209](https://doi.org/10.1371/journal.pcbi.1010209)
