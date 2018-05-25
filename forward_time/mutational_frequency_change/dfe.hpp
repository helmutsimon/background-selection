#ifndef DFE_HPP
#define DFE_HPP

#include <vector>
#include <cmath>
using namespace std;

#include "stats.hpp"
#include "mutation.hpp"


class DeltaDFE{
public:
	double w;
	double U;

	DeltaDFE(double s, double U): w{std::exp(s)}, U(U) {};

	double get_mutation_rate() { return U; }

	double get_fitness_effect(Random & random){ return w; }

	Mutation get_mutation(Random & random, LabelGenerator & label_generator){
		return Mutation{label_generator.get_next_label(), this->get_fitness_effect(random)};
	}

	void mutate_individual(Random & random, LabelGenerator & label_generator, Individual & individual){
		int no_mutations = sample_poisson(random,U);
		if (no_mutations > 0){
			for (int i = 0; i < no_mutations; ++i){
				individual.add_mutation(get_mutation(random,label_generator));
			}
		}
	}

	void mutate_individual(Random & random, Individual & individual){
		int no_mutations = sample_poisson(random,U);
		if (no_mutations > 0){
			for (int i = 0; i < no_mutations; ++i){
				individual.fitness *= get_fitness_effect(random);
			}
		}
	}

};




#endif