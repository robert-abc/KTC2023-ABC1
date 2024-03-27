# Post-processing electrical impedance tomography reconstructions with incomplete data using convolutional neural networks (Source code)

## Brief description of our algorithm in the context of Kuopio Tomography Challenge 2023 (KTC2023) [[1]](#1)

We based our algorithm the smoothness prior provided by the KTC2023 organising committee [[2]](#2). We modified it to include post-processing steps, including one with a convolutional neural network (CNN), to remove artifacts from electrode disconnection and to segment the results. The official KTC2023 results are available [here](https://www.fips.fi/KTCresults.php)

Here, we include two codes: The original proposal submitted to the challenge (Matlab) and the proposal using hyperparameter optimization (Python). Each subfolder includes installation and usage instructions. We also included the weights of the CNN after the training.

We used Google Colaboratory (Colab) to train the CNN. The training set and the Python code for the network architecture and optimization are available at Zenodo [here](https://zenodo.org/records/10801591) (DOI: 10.5281/zenodo.10801591)

Please, contact us if you need further information. 

## Authors
* Roberto Gutierrez Beraldo¹* - roberto.gutierrez@ufabc.edu.br
* Leonardo Alves Ferreira¹ - leonardo.alves@ufabc.edu.br
* Fernando Silva de Moura² - fernando.moura@ufabc.edu.br
* André Kazuo Takahata¹ - andre.t@ufabc.edu.br
* Ricardo Suyama¹ - ricardo.suyama@ufabc.edu.br
  
*Corresponding author

¹Federal University of ABC - [UFABC](https://www.ufabc.edu.br/) - Campus Santo André - Avenida dos Estados, 5001 - Bairro Bangu, Santo André - CEP: 09280-560 (Brazil)

²Federal University of ABC - [UFABC](https://www.ufabc.edu.br/) - Campus São Bernardo - Alameda da Universidade, s/nº - Bairro Anchieta, São Bernardo do Campo - CEP: 09606-405 (Brazil)

## References

<a id="1">[1]</a> 
M. Räsänen et al.
“Kuopio Tomography Challenge 2023”. [Online]. Available at: https://www.fips.fi/KTC2023.php.

<a id="2">[2]</a> 
M. Räsänen et al. 
“Kuopio Tomography Challenge 2023 open electrical impedance tomographic dataset (KTC 2023)”. Available at: https://zenodo.org/records/8252370.
