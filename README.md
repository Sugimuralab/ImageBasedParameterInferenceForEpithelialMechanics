Parameter estimation and model selection for epithelial mechanics
===

## Description

This is a script for parameter esitmation and model selection for epithelial mechanics, proposed in Ogita et al. 2022 [1]. 

The most predictive tension model of cell jucntions and its model parameter will be obtained from the cell geometry data of epithelial tissue.
Five models (model A-E) will be considered.

See Ogita et al. 2022[1] for details.

The script for Bayesian force inference is [here (https://github.com/IshiharaLab/BayesianForceInference).

## Requirement

* statsmodels 0.12.2


## Usage

1. Prepare an input file in the same format as the attached sample (Sample/sample/).
2. Change the variable "filename" in ParameterEstimation.py to the input file in step 1.
3. Run ParameterEstimation.py on IDE or IPython.

## Licence

[MIT](https://github.com/tcnksm/tool/blob/master/LICENCE)

## Reference

1. Goshi Ogita, Takefumi Kondo, Keisuke Ikawa, Tadashi Uemura, Shuji Ishihara & Kaoru Sugimura (2022)<br>"Image based parameter inference for epithelial mechanics"
