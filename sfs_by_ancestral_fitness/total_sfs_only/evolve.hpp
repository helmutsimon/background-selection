#ifndef EVOLVE_HPP
#define EVOLVE_HPP
#include <cmath>

#include "stats.hpp"
#include "individual.hpp"
#include "dfe.hpp"
#include "diversity.hpp"


class Measurement {
    public:
        FrequencySpectrum frequency_spectrum;
};

typedef std::vector<Measurement> Results;

template<class NonNeutralDFE, class NeutralDFE>
inline Results evolve(Random & random, LabelGenerator & label_generator, int N, 
						NonNeutralDFE & non_neutral_dfe, NeutralDFE & neutral_dfe, int total_samples)
{
	// declare stuff we want to track in this forward time simulation
	Results results;

	MutationList fixed_mutations;
	FrequencySpectrum frequency_spectrum(0.0,N - 1);

	// initialize populations
	//need two populations to avoid memory reallocation
	Population population; // this timestep
	Population new_population; //next timestep
	population.reserve(N);
	new_population.reserve(N);

	double total_fitness = 0.0;
	for (int i = 0; i < N; i++){
	    population.push_back(Individual{total_fitness+=1.0});
	}

	int t(0);
	int collected_samples(0);
	bool at_equilibrium = false;
	int output_interval(0);
	while (collected_samples < total_samples)
	{	
		// resample population
		double new_total_fitness = 0.0;
		for (int i = 0; i < N; i++){
			new_population.push_back(*std::lower_bound(population.begin(), population.end(), sample_uniform(random)*total_fitness));
			
			non_neutral_dfe.mutate_individual(random,new_population.back());
			neutral_dfe.mutate_individual(random,label_generator,new_population.back());
			
			new_total_fitness += new_population.back().fitness;
			new_population.back().weight = new_total_fitness;
		}

		population.swap(new_population);
		total_fitness = new_total_fitness;
		new_population.clear();

		// look for fixed mutations and record them
		bool found_fixed = true;

		while (found_fixed){
			auto query_for_fixation = population.front().get_earliest_mutation();
			// std::cout << "Querying mutation " << query_for_fixation.label << std::endl;

			if (query_for_fixation.label != null_label){
				// check whether all lineages have this mutation, if not it's not fixed
				for (auto & individual: population){
					if (individual.get_earliest_mutation().label != query_for_fixation.label){
						found_fixed = false;
						break;
					}
				}	
				// if it's fixed, remove it from tracked mutations in population and put in fixed_mutations;
				if (found_fixed){
					fixed_mutations.push_back(query_for_fixation);
					for (auto & individual: population){
						if (individual.get_earliest_mutation().label == query_for_fixation.label)
							individual.remove_earliest_mutation();
					}
				}		
			}
			else 
				found_fixed = false; 
		}

		if (at_equilibrium) {
			if (t % output_interval == 0){
					collected_samples++;

					Population sample;
					sample.reserve(N);
					for (auto & individual:population){
						sample.push_back(individual);
					}
					
					frequency_spectrum += calculate_frequency_spectrum(sample);

					FrequencySpectrum measurement(frequency_spectrum);
					for (auto & frequency:measurement)
					{
						frequency /= collected_samples;
					}

					results.push_back(Measurement{measurement});
			}
		}	

		if ((fixed_mutations.size()>0) && !(at_equilibrium)){ 
			output_interval = t;
			at_equilibrium = true;
		}

		t++;
	}
	return results;
}


#endif
