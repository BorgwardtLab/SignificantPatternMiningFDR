# SignificantPatternMiningFDR

Code and Data for the paper: *FASM and FAST-YB: Significant Pattern Mining with False Discovery Rate Control* (ICDM 2023).

In significant pattern mining, i.e. the task of dis- covering structures in data that exhibit a statistically significant association with class labels, it is often needed to have guarantees on the number of patterns that are erroneously deemed as statistically significant by the testing procedure. A desirable property, whose study in the context of pattern mining has been limited, is to control the expected proportion of false positives, often called the false discovery rate (FDR). 
This repository contains the code for two novel algorithms for mining statistically significant patterns under FDR control. 

The first one, **FASM**, builds upon the Benjamini-Yekutieli procedure and exploits the discrete nature of the test statistics to increase its computational efficiency and statistical power. 
The second one, **FAST-YB**, incorporates the Yekutieli-Benjamini permutation testing procedure to account for interdependencies among patterns, which allows for a further increase in statistical power. 

### Citing our work
FASM and Fast-YB are described in the following paper:
> Paolo Pellizzoni and Karsten Borgwardt. *FASM and FAST-YB: Significant Pattern Mining with False Discovery Rate Control*. ICDM 2023. 

### Appendix to the main paper
The file ```Appendix.pdf``` contains details, proofs and experimental results that were omitted from the main paper due to space constraints.

### Organization of the codebase
The folders ```fasm/``` and ```fastYB/``` contain respectively the source codes for FASM and Fast-YB.
The folder ```wy_light/``` contains the original implementations of LAMP and WY-Light.
All source codes can be compiled with make and gcc.

The datases we used are contained in the ```datasets.zip``` archive, which is obtained from https://github.com/VandinLab/TopKWY.
