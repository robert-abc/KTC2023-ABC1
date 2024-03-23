# Proposed method installation and usage instructions (MATLAB)

## Method installation and requirements
* The Matlab files available in this subfolder.
* Most files are similar or identical to the original KTC codes. The differences are:
    * We added the CNN file ("ultimate_cnn1.h5")
    * We modified main.m to include the post-processing steps

| Package | Version |
| ------------- | ------------- |
| Matlab | R2023a Update 5 | 
| Deep Learning Toolbox | 14.6 | 
| Image Processing Toolbox | 11.7 | 

Regarding the external codes and toolboxes: 
* We used the function importKerasLayers from the Deep Learning Toolbox Converter for TensorFlow Models (V.23.1.0).
    * It is necessary to download and install it from https://www.mathworks.com/matlabcentral/fileexchange/64649-deep-learning-toolbox-converter-for-tensorflow-models
    * More info at: https://www.mathworks.com/help/deeplearning/ref/importkeraslayers.html
    * When importing the CNN, Matlab shows the following warning, but imports the CNN nevertheless:
        * Warning: File 'ultimate_cnn1.h5' was saved in Keras version '2.11.0'. Import of Keras versions newer than '2.6.0' has not been thoroughly tested yet. The imported model might not exactly match the model saved in the Keras file.
    * Note that a new function importNetworkFromTensorFlow has been available since Matlab R2023b, but we did not use it.
* We based our implementation of the soft thresholding on the following function from the UNLocBoX toolbox:
    * https://epfl-lts2.github.io/unlocbox-html/doc/utils/soft_threshold_code.html
    * In this case, there's no need to download this toolbox or this function


## Usage instructions and example: Running with a callable function from the command line

By the rules of the KTC2023, it was expected a main routine with three arguments: 
* *Your main routine must require three input arguments:*
1. *(string) Folder where the input files are located (.mat files)*
1. *(string) Folder where the output files must be stored*
1. *(int) (int) Group category number (Values between 1 and 7)*

To run the code:
* Install Deep Learning Toolbox Converter for TensorFlow Models
* Download the files from the folder \KTC2023_Codes_1_Github
* Create the folder where the input files are located
* Create the folder where the output files must be stored
* Run the Proposal_function.m, which calls the main.m, using the three arguments from above appropriately.
* Run the scoring function
    * Although the scoringFunction.m is included in the folder, we do not evaluate the scores for each output. This part of the code is given as comments in the main.m function.
* Note that, inside main.m, there is the parameter PAUSE_FOR_EACH_TARGET:
  * If 1, we plot images from all steps (default).
  * If 0, there is no plot. 
