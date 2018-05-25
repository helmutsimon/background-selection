Forward-time Wright-Fisher simulations are implemented using an approach that is closely modeled after that of Good et. al. PLoS Genetics (2014).

  - **sfs_by_ancestral_fitness/**: scripts in this directory take a set of population genetic parameters as the input, and output the total site frequency spectrum of the population, as well as marginal site frequency spectra of mutations founded in individuals of different fitnesses. Scripts in the sub-directory **sfs_by_ancestral_fitness/total_sfs_only/** output only the total site frequency spectrum.
   - **mutational_frequency_change/**: scripts in this directory output pairs of frequencies of all mutations in the population at two different timepoints.

### Requirements and Installation:
The only requirement is a C++ compiler that supports the C++11 standard. To install, simply compile the "simulate.cpp" source file and save it as 'simulate_frequency_spectrum' or 'simulate_two_timepoint_sampling', depending on the type of simulation.

### Usage:

Scripts take a number of required command line arguments specifying the population genetic parameters, as well as the number of samples to average over. The usage for the two types of simulations is: 

	./simulate_frequency_spectrum N NUd NUn Ns num_samples

	./simulate_two_timepoint_sampling N NUd NUn Ns time_interval num_samples

where 
- ``N`` is the size of the population, 
- ``NUd`` and ``NUn`` specify the product of the population size and the deleterious and neutral mutation rates, respectively, 
- ``Ns`` specifies the product of the population size and the selective cost of deleterious mutations, and 
- ``time_interval`` represents the time in generations between the two samples.
