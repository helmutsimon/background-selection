#ifndef EVOLVE_HPP
#define EVOLVE_HPP
#include <cmath>
#include <set>

#include "stats.hpp"
#include "individual.hpp"
#include "dfe.hpp"


class Results {
    public:
        // int time;        
        // double average_fitness;
        std::vector<std::pair<int,int>> counts;
};

template<class NonNeutralDFE, class NeutralDFE>
inline Results evolve(Random & random, LabelGenerator & label_generator, int N, 
						NonNeutralDFE & non_neutral_dfe, NeutralDFE & neutral_dfe, int total_samples, int time_interval)
{
	// declare stuff we want to track in this forward time simulation
	Results results;

	std::set<Mutation> fixed_mutation_set;
	MutationList mutations;
	std::set<Mutation> seen_mutations;
	std::set<Mutation> unique_mutations;
	std::multiset<Mutation> all_mutations;

	// interval between which mutations are tracked
	// int time_interval = 2 * static_cast<int>((log(non_neutral_dfe.U/(-log(non_neutral_dfe.w))))/(-log(non_neutral_dfe.w)));

	// std::cerr << "time interval is " << time_interval << std::endl;
	// initialize populations
	Population population; // this timestep
	Population new_population; //next timestep
	Population sample;
	population.reserve(N);
	new_population.reserve(N);
	sample.reserve(N);

	double total_fitness = 0.0;
	for (int i = 0; i < N; i++){
	    population.push_back(Individual{total_fitness+=1.0});
	}

	int t(0);
	int collected_samples(0);
	bool at_equilibrium = false;
	int output_interval(0);
	double max_fitness = 1.0;
	while (collected_samples < total_samples)
	{	
		// resample population
		double new_total_fitness = 0.0;
		double new_max_fitness = 0.0;
		for (int i = 0; i < N; i++){
			new_population.push_back(*std::lower_bound(population.begin(), population.end(), sample_uniform(random)*total_fitness));
			
			non_neutral_dfe.mutate_individual(random,new_population.back());
			neutral_dfe.mutate_individual(random,label_generator,new_population.back());
			
			if (new_population.back().fitness > new_max_fitness)
				new_max_fitness = new_population.back().fitness;

			new_total_fitness += new_population.back().fitness;
			new_population.back().weight = new_total_fitness;
		}
		
		population.swap(new_population);
		total_fitness = new_total_fitness;
		max_fitness = new_max_fitness;
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
					fixed_mutation_set.insert(query_for_fixation);
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
					seen_mutations.clear(); mutations.clear();
					all_mutations.clear(); unique_mutations.clear();
					// compile all mutations and their frequencies
					for(auto & individual: population){
						mutations.insert(mutations.end(),individual.mutations->begin(), individual.mutations->end());
					}
					unique_mutations.insert(mutations.begin(), mutations.end());
					all_mutations.insert(mutations.begin(), mutations.end());
					for (auto & mutation : unique_mutations){
						seen_mutations.emplace_hint(seen_mutations.end(),mutation.label,0.0,all_mutations.count(mutation));
					}
			}

			if ((t > 2*output_interval) && (t% output_interval == time_interval)){
				collected_samples++;
				// std::cerr << "collected sample " << collected_samples << std::endl;
				// compile mutations in the sample
				mutations.clear(); all_mutations.clear();

				for(auto & individual: population){
						mutations.insert(mutations.end(),individual.mutations->begin(), individual.mutations->end());
				}
				all_mutations.insert(mutations.begin(), mutations.end());
				for (auto & mutation: seen_mutations){
					int counts_now = all_mutations.count(mutation);
					
					if (counts_now == 0){
						if (fixed_mutation_set.count(mutation))
							counts_now = N;
					}
					std::pair<int,int> new_pair(mutation.f_initial,counts_now);
					results.counts.push_back(new_pair);
				}

			}
		}	

		if ((fixed_mutation_set.size()>0) && !(at_equilibrium)){ 
			output_interval = t;
			// std::cerr << "output_interval is " << output_interval << std::endl;
			if (output_interval < time_interval)
				output_interval = 2*time_interval;
			at_equilibrium = true;
		}

		t++;
	}
	return results;
}


#endif