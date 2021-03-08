# Maize Machine Learning Models
In this project, we assess the performance of three genomic prediction methods: BayesA, rrBLUP, and neural network. We predict three traits: anthesis GDD, kernel row number, and plant height.

# Project Objectives:
Predict three maize phenotypic traits with genomic markers using three different methods.
Quantify the performance of each method and assess their strengths and weaknesses.
Develop clean, well documented Python code.
Use GitHub to facilitate collaboration.

# Dependencies
The following python libraries are required.
* `numpy`
* `pandas`
* `seaborn`

# Installation and Usage:
For the BayesA and rrBLUP files, several file paths are hard-coded into the scripts and you will have to modify these.

To run the neural network code, run the script as follows.
To run the code: `python nn_pred.py US_NAM`
To display a help file: `python nn_pred.py -h`

To generate graphs make sure you `cd` into the *src* folder. Then run `python graphs_nn_pred.py Results/`
