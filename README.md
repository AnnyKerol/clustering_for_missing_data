# Adaptive kernel fuzzy clustering for missing data

In this work, the Fuzzy C-means clustering algorithm is used considering the kernelization of the metric with local adaptive distances (VKFCM-K-LP) under three types of strategies to deal with missing data, Whole Data Strategy (WDS), Partial Distance Strategy (PDS) and Optimal Completion Strategy (OCS).

#### Step by step to execute VKFCM-K-LP with  WDS, PDS and OCS strategies
Example with the iris dataset with 5% missings
1. In a file with a .sh extension, the following arguments are passed <br>
 **./vkfcm-k-lp irisNA5.data r-iris5.txt iris.par iris.sig idx-iris5.txt**
   - vkfcm-k-lp: run the executable on linux for example
   - irisNA5.data: missing dataset
   - r-iris5.txt: is the name of the file that will store the results
   - iris.par: is the simulation configuration file
   - iris.sig: is a file with sigma values
   - idx-iris5.txt: is a file to store the values of the CR, F-Measure and OERC indexes of each repetition.
2. The iris.par file must have the following structure
   **150 5 3 2.0 300 1e-10 100**
   - 150: is the number of observations
   - 5: is the number of variables
   - 2.0: is the fuzzy exponent
   - 300: is the maximum number of iterations of the algorithm
   - 1e-10: is the tolerance for convergence
   - 100: is the number of repetitions.
3. To execute the OCS strategy it is necessary to add the **bdimput5.txt** argument in the .sh file <br>
   **./vkfcm-k-lp irisNA5.data r-iris5.txt iris.par iris.sig idx-iris5.txt bdimput5.txt**
   - bdimput5.txt: is the complete dataset with missings imputed by the OCS strategy.
