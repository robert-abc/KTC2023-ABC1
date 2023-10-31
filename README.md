# Post-processing smoothness prior EIT reconstruction using CNNs

## Brief description of our algorithm in the context of Kuopio Tomography Challenge 2023 (KTC2023) [[1]](#1)

Our algorithm is based on the smoothness prior provided by the KTC2023 organising committee [[2]](#2). We modified it to include post-processing steps, including one with a CNN, to extract meaninful information and to segment the resuts. This is our first submission to the KTC2023.

A brief and visual overview of our proposal is attached in the [Visual read me PDF](visual_readme.pdf). It also includes examples of reconstructions using data from the traning set. 

## Authors
* Roberto Gutierrez Beraldo¹* - roberto.gutierrez@ufabc.edu.br
* Leonardo Ferreira Alves¹ - leonardo.alves@ufabc.edu.br
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
    The one line code that calls main.m is the Proposal_function.m
* We trained the CNN using Tensorflow/keras in Google Colaboratory (Colab). We did not include the python file, the code to generate the traning pairs, or the training pairs themselves here. Please, contact us if you need further information other the [Visual read me](visual_readme.pdf).


| Package | Version |
| ------------- | ------------- |
| Matlab | R2023a Update 5 | 
| Deep Learning Toolbox | 14.6 | 
| Image Processing Toolbox | 11.7 | 

### External codes

* We used the function importKerasLayers from the Deep Learning Toolbox Converter for TensorFlow Models (V.23.1.0).
    * It is necessary to download and install it from https://www.mathworks.com/matlabcentral/fileexchange/64649-deep-learning-toolbox-converter-for-tensorflow-models
    * More info at: https://www.mathworks.com/help/deeplearning/ref/importkeraslayers.html
    * Note that a new function importNetworkFromTensorFlow is available since Matlab R2023b, but we did not use it.
* We based our implementation of the soft thresholding in the following function from the UNLocBoX toolbox:
    * https://epfl-lts2.github.io/unlocbox-html/doc/utils/soft_threshold_code.html
    * In this case, there's no need to download this toolbox or this function


### Usage instructions and example: Running with a callable function from the command line

By the rules of the KTC2022, it was expected a one-line command: 
* *Your main routine must require three input arguments:*
1. *(string) Folder where the input image files are located*
1. *(string) Folder where the output images must be stored*
1. *(int) Difficulty category number. Values between 1 and 7*

After downloading the files from the folder \KTC2023_Codes_1_Github, ...

## References

<a id="1">[1]</a> 
M. Räsänen et al.
“Kuopio Tomography Challenge 2023”. [Online]. Available at: https://www.fips.fi/KTC2023.php.

<a id="2">[2]</a> 
M. Räsänen et al. 
“Kuopio Tomography Challenge 2023 open electrical impedance tomographic dataset (KTC 2023)”. Available at: https://zenodo.org/records/8252370.
