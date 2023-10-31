# Post-processing smoothness prior EIT reconstruction using CNNs

## Brief description of our algorithm in the context of Kuopio Tomography Challenge 2023 (KTC2023) [[1]](#1)

Our algorithm is based on the smoothness prior provided by the KTC2023 organising committee [[2]](#2). We modified it to include post-processing steps, including one with a CNN, to extract meaninful information and to segment the resuts. This is our first submission to the KTC2023.

A brief and visual overview of our proposal is attached in the pdf: [A visual read me](visual_readme.pdf). It also includes examples of reconstructions using data from the traning set. 

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
* The Python codes are available in this repository, see main.py and the /utils folder.

* We ran our codes using Google Colaboratory (Colab), but it results in a large list of packages (obtained by pip freeze > requirements.txt) and not all of them are necessary.
* It is possible to create an Anaconda environment "by hand" given the packages list. In the following table, there is a small list of the main packages we used (with "import").

| Package | Version |
| ------------- | ------------- |
| Python | 3.10.11 | 
| Numpy | 1.22.4 | 
| Matplotlib | 3.7.1 | 
| Scipy | 1.10.1 | 
| Skimage | 0.19.3 |
| Pillow | 8.4.0 | 
| Torch | 2.0.0+cu118 | 
| ODL | 1.0.0 | 

### Usage instructions and example: Running with a callable function from the command line

By the rules of the HTC 2022, it was expected a one-line command: 
* *Your main routine must require three input arguments:*
1. *(string) Folder where the input image files are located*
1. *(string) Folder where the output images must be stored*
1. *(int) Difficulty category number. Values between 1 and 7*
* *Python: The main function must be a callable function from the command line. 

After the setup, it is possible to run our code following these rules. Considering the difficulty group 7: 
* python main.py 'example/input' 'example/output' 7

See, for instance, the Section "Generating results" from the example notebook [Here](/notebook_example.ipynb).


## References

<a id="1">[1]</a> 
M. Räsänen et al.
“Kuopio Tomography Challenge 2023”. [Online]. Available at: https://www.fips.fi/KTC2023.php.

<a id="2">[2]</a> 
M. Räsänen et al. 
“Kuopio Tomography Challenge 2023 open electrical impedance tomographic dataset (KTC 2023)”. Available at: https://zenodo.org/records/8252370.
