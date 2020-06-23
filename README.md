# A reduced universum twin support vector machine for class imbalance learning (RUTSVM-CIL)
This is implementation of the paper: B. Richhariya, M. Tanveer, A reduced universum twin support vector machine for class imbalance learning, Pattern Recognition, Volume 102, 2020, 107150, https://doi.org/10.1016/j.patcog.2019.107150.

Description of files:

main.m: selecting parameters of RUTSVM-CIL using k fold cross-validation. One can select parameters c, mu and e to be used in grid-search method.

rutsvm.m: implementation of RUTSVM-CIL algorithm. Takes parameters c, mu, e, and training data and test data, and provides accuracy obtained and running time.

For quickly reproducing the results of the RUTSVM-CIL algorithm, we have made a newd folder containing a sample dataset. One can simply run the main.m file to check the obtained results on this sample dataset. To run experiments on more datasets, simply add datasets in the newd folder and run main.m file. The code has been tested on Windows 10 with MATLAB R2017a.

This code is for non-commercial and academic use only.
Please cite the following paper if you are using this code.

Reference: B. Richhariya, M. Tanveer, A reduced universum twin support vector machine for class imbalance learning, Pattern Recognition, Volume 102, 2020, 107150, https://doi.org/10.1016/j.patcog.2019.107150.

For further details regarding working of algorithm, please refer to the paper.

If further information is required you may contact on: phd1701241001@iiti.ac.in.
