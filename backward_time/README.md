These implementations of the structured coalescent and the Bolthausen-Sznitman coalescent are provided only for completeness and are reproduced with minor modifications from the scripts accompanying Good et. al. _PLoS Genetics_ (2014).

### Requirements and Installation:
The only requirement is a C++ compiler that supports the C++11 standard.

### Usage:

Scripts take a number of required command line arguments specifying the population genetic parameters, as well as the number of samples to average over. The usage for the two types of simulations is: 

	./struct NUd Ns n num_samples

	./bsc n num_samples

where 
- ``NUd`` specify the product of the population size and the deleterious mutation rates, 
- ``Ns`` specifies the product of the population size and the selective cost of deleterious mutations,  
- ``n`` specifies the sample size, and
- ``num_samples`` specifies the number of genealogies to sample.