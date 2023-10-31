# Post-processing smoothness prior EIT reconstruction using CNNs

## Brief description of our algorithm in the context of Kuopio Tomography Challenge 2023 (KTC2023) [[1]](#1)

Our algorithm is based on the smoothness prior provided by the KTC2023 organising committee [[2]](#2). We modified it to include post-processing steps, including one with a CNN, to extract meaninful information and to segment the resuts. This is our first submission to the KTC2023.

A brief and visual overview of our proposal is attached in the pdf: [A visual read me](visual_readme.pdf). It also includes examples of reconstructions from the traning set. 

## Authors
* Roberto Gutierrez Beraldo¹* - roberto.gutierrez@ufabc.edu.br
* Leonardo Ferreira Alves¹ - leonardo.alves@ufabc.edu.br
* Fernando Silva de Moura¹ - fernando.moura@ufabc.edu.br
* André Kazuo Takahata¹ - andre.t@ufabc.edu.br
* Ricardo Suyama¹ - ricardo.suyama@ufabc.edu.br
* 
*Corresponding author

¹Federal University of ABC - https://www.ufabc.edu.br/.

-Avenida dos Estados, 5001 - Bairro Bangu, Santo André - CEP: 09280-560 (Brazil)

-Alameda da Universidade, s/nº - Bairro Anchieta, São Bernardo do Campo - CEP: 09606-405 (Brazil)



## 1. HTC 2022 [[3]](#3): Overview, rules, and constraints
Challenge URL: https://www.fips.fi/HTC2022.php

* The challenge consists in reconstructing limited-angle computed tomography (CT) images,  i.e. reducing the number of the total projections to a region,  although it is expected to be a general-purpose (CT reconstruction) algorithm.  
* Several design choices of the algorithm were related to the challenge. 
* After the tomographic image is reconstructed, it is necessary to segment it into two parts, binarizing it: the acrylic disk (1's) and the background, including holes (0's).   
* For more details, please refer to https://www.fips.fi/Helsinki_Tomography_Challenge_2022_v1.pdf.
* The training dataset consists of only 5 sinograms, while the test set contains 21 cases divided in 7 levels of difficulty. They are avaliable at: https://zenodo.org/record/7418878. 
* We know the subsampling of each difficulty group, that is: The information given includes the initial (random) angle, the angular range, and the angular increment. The sinogram size varies for each difficulty group. 
* We do not have information for each difficulty group about the number and size of the holes. 

## 2. Forward problem 
We consider the forward problem, i.e. calculating the sinogram from the tomographic image, as 

$$Ax = y$$

where $x$ is the tomographic image, $A$ is the (linear) forward model, and $y$ is the resulting sinogram [[4]](#4).  

* The cone beam computed tomography considers the detector is flat. As the image is 2D, the cone beam CT can be approximated by the fan beam CT. 
* We use the ODL Pytorch extension (https://github.com/odlgroup/odl) to obtain $A$ not as an explicit matrix, but as an object.
* By the instructions, the submitted algorithm does not need to subsample the test data (as this has already been done). In this way, we only need to create an appropriate $A$ for each difficulty group, given the initial angle and the angular range. 


### Note:
* The above $A$ matrix size is calculated considering an image of 512x512 = 262144 values and  considering 560 (detector points) x [360*2+1] (a projection for each 0.5° angular increment) = 403760 values when including all angles available. Although sparse, $A$ has 262144 x 403760 values considering full-angle CT, so it is often not practical to explicitly calculate $A$.  
 


## 3. Proposed method installation, usage instructions, and examples

### 3.1 Method installation and requirements
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

### 3.2 Usage instructions and example: Running with a callable function from the command line

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
