# Post-processing smoothness prior EIT reconstruction using CNNs

## Brief description of our algorithm in the context of Kuopio Tomography Challenge 2023 (KTC2023) [[1]](#1)

Our algorithm is based on the smoothness prior provided by the KTC2023 organising committee [[2]](#2). We modified it to include post-processing steps, including one with a CNN, to extract meaningful information and to segment the results. This is our first submission to the KTC2023.

A brief and visual overview of our proposal is attached in the [Visual read me PDF](visual_readme.pdf). It also includes examples of reconstructions using data from the training set. 

## Authors
* Roberto Gutierrez Beraldo¹* - roberto.gutierrez@ufabc.edu.br
* Leonardo Alves Ferreira¹ - leonardo.alves@ufabc.edu.br
* Fernando Silva de Moura² - fernando.moura@ufabc.edu.br
* André Kazuo Takahata¹ - andre.t@ufabc.edu.br
* Ricardo Suyama¹ - ricardo.suyama@ufabc.edu.br
  
*Corresponding author

¹Federal University of ABC - [UFABC](https://www.ufabc.edu.br/) - Campus Santo André - Avenida dos Estados, 5001 - Bairro Bangu, Santo André - CEP: 09280-560 (Brazil)

²Federal University of ABC - [UFABC](https://www.ufabc.edu.br/) - Campus São Bernardo - Alameda da Universidade, s/nº - Bairro Anchieta, São Bernardo do Campo - CEP: 09606-405 (Brazil)

## Proposed method installation and usage instructions

### Method installation and requirements
* The Matlab files available in this repository, see /KTC2023_Codes_1_Github folder.
* Most files are similar or identical to the original KTC codes. The differences are:
    * We added the CNN file ("ultimate_cnn1.h5")
    * We modified main.m to include the post-processing steps
* We trained the CNN using Tensorflow/keras in Google Colaboratory (Colab). We did not include the Python file, the code to generate the training pairs or the training pairs themselves here. Please, contact us if you need further information.

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


### Usage instructions and example: Running with a callable function from the command line

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

## References

<a id="1">[1]</a> 
M. Räsänen et al.
“Kuopio Tomography Challenge 2023”. [Online]. Available at: https://www.fips.fi/KTC2023.php.

<a id="2">[2]</a> 
M. Räsänen et al. 
“Kuopio Tomography Challenge 2023 open electrical impedance tomographic dataset (KTC 2023)”. Available at: https://zenodo.org/records/8252370.
