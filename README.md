# rEMT
+ In this repository, we provide Python and R codes for reproducing the analyses of Boolean network model shown in the following reference paper.  
+ The code was originally written by Namhee Kim and has been developed with the help of many others (Jonghoon Lee, etc.).  

**Reference paper:**  A Cell-Fate Reprogramming Strategy Reverses Epithelial-to-Mesenchymal Transition of Lung Cancer Cells While Avoiding Hybrid States
https://doi.org/10.1158/0008-5472.CAN-22-1559

## Requirements

+ The code is developed using Python 3.7.x and R 4.2.x
+ Some packages must be installed:
  
  pyboolnet: <https://github.com/hklarner/pyboolnet>  
  BoolNet: 10.32614/CRAN.package.BoolNet


## Implementation

+ Import networks in BoolNet format (from the `network` directory)  
+ A simple toy network example is provided.  

+ 1-Attractor landscape analysis.ipynb : This is an example of attractor simulation & analysis of Boolean network model.  
+ 2-Molecular state ambiguity.ipynb : This is an example of how Boolean network model was analyzed based on the concept of frustration.  
+ 3~8-Jupyter Notebook files : Note that these files reproduce our results shown in the main figures.  


## Notes

+ Main functions are in the `modules` directory.  
+ All outputs (e.g., simulation results) are already provided in the `result` directory.  
+ Some results can be obtained after adjusting certain variables in the script.  
