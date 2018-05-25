#ifndef DIVERSITY_HPP
#define DIVERSITY_HPP

#include <valarray>
#include <set>
#include <set>
#include "mutation.hpp"
#include "individual.hpp"

typedef std::valarray<double> FrequencySpectrum;

FrequencySpectrum calculate_frequency_spectrum(Population & sample){
	int n = sample.size();
	FrequencySpectrum frequency_spectrum(0.0,n-1);

	MutationList mutations;
	for(auto & individual: sample){
		mutations.insert(mutations.end(),individual.mutations->begin(), individual.mutations->end());
	}

 	std::set<Mutation> unique_mutations(mutations.begin(), mutations.end());
 	std::multiset<Mutation> all_mutations(mutations.begin(), mutations.end());

	for (auto mutation : unique_mutations){
		int i = all_mutations.count(mutation);
		if (i < n)
			frequency_spectrum[i-1]++;
	}

	return frequency_spectrum;
}

#endif