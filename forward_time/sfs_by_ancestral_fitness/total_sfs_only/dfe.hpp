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

template <class DFE1, class DFE2>
class CompositeDFE{
public:
	DFE1 first_dfe;
	DFE2 second_dfe;

	CompositeDFE(DFE1 first_dfe, DFE2 second_dfe): first_dfe(first_dfe), second_dfe(second_dfe) {};

	double get_mutation_rate() { return first_dfe.get_mutation_rate() + second_dfe.get_mutation_rate(); }

	double get_fitness_effect(Random & random){
		double p = sample_uniform(random);
		
		if (p < first_dfe.get_mutation_rate()/(first_dfe.get_mutation_rate() + second_dfe.get_mutation_rate()))
			return first_dfe.get_fitness_effect(random);
		else
			return second_dfe.get_fitness_effect(random);
	};

	Mutation get_mutation(Random & random, LabelGenerator & label_generator){
		Mutation new_mutation{label_generator.get_next_label(), get_fitness_effect(random)};
		return new_mutation;
	}

	void mutate_individual(Random & random, LabelGenerator & label_generator, Individual & individual){
		first_dfe.mutate_individual(random,label_generator,individual);
		second_dfe.mutate_individual(random,label_generator,individual);
	}

	void mutate_individual(Random & random, Individual & individual){
		first_dfe.mutate_individual(random,individual);
		second_dfe.mutate_individual(random,individual);
	}

};



#endif