#ifndef DIVERSITY_HPP
#define DIVERSITY_HPP

#include <valarray>
#include <set>
#include <set>
#include "mutation.hpp"
#include "individual.hpp"


inline bool mut_comp(const Mutation & x, const Mutation & y){return x.ancestral_background_fitness < y.ancestral_background_fitness;}
inline bool operator < (const Mutation & x, const int y){return x.ancestral_background_fitness < y;}
inline bool operator < (const int x, const Mutation & y){return x < y.ancestral_background_fitness;}

typedef std::valarray<double> FrequencySpectrum;

void divide(std::valarray<double> v_array, double x){
	for (auto item:v_array){
		item /= x;
	}
}

FrequencySpectrum calculate_frequency_spectrum(MutationList & mutations, int n){
	FrequencySpectrum frequency_spectrum(0.0,n-1);
	std::set<Mutation> unique_mutations(mutations.begin(), mutations.end());
 	std::multiset<Mutation> all_mutations(mutations.begin(), mutations.end());
	for (auto mutation : unique_mutations){
		int i = all_mutations.count(mutation);
		if (i < n){
			frequency_spectrum[i-1]++;
		}
	}
	return frequency_spectrum;
}

FrequencySpectrum calculate_frequency_spectrum(Population & sample){
	int n = sample.size();

	MutationList mutations;
	for(auto & individual: sample){
		mutations.insert(mutations.end(),individual.mutations->begin(), individual.mutations->end());
	}

	return calculate_frequency_spectrum(mutations, n);
}

std::vector<FrequencySpectrum> calculate_frequency_spectra_for_each_background(Population & sample, int no_categories, double w){
	
	std::vector<FrequencySpectrum> results;
	results.reserve(no_categories);
	int n = sample.size();

	// collect all mutations in sample
	MutationList mutations;
	for(auto & individual: sample){
		mutations.insert(mutations.end(),individual.mutations->begin(), individual.mutations->end());
	}

	int num_mutants;

	// sort on background fitness: mutations will be sorted in ascending order
	std::sort(mutations.begin(),mutations.end(),mut_comp);

	num_mutants = mutations.size();
	// std::cout << "total\t" << num_mutants << std::endl;
	// for (auto & mut:mutations){
	// 	std::cout << mut.ancestral_background_fitness << " ";
	// }
	// std::cout << std::endl;
	// prepare a container for subsets of mutations
	MutationList subset;
	subset.reserve(mutations.size());


	int num_accounted = 0;

	// calculate the conditional site frequency spectra for mutations that originated on backgrounds with < no_categories-1 mutations
	for (int i = 0; i<no_categories-1; i++){
		subset.clear();
		int target_w = -i;

		auto lower = std::lower_bound(mutations.begin(), mutations.end(), target_w);
   		auto upper = std::upper_bound(mutations.begin(), mutations.end(), target_w);
   		num_accounted += upper - lower;
   		
   		if (upper > lower){
			subset.insert(subset.end(),lower, upper);
			results.push_back(calculate_frequency_spectrum(subset, n));
		}
		else
			results.push_back(FrequencySpectrum(0.0,n-1));
	}
	// now calculate the frequency spectrum
	subset.clear();
	auto upper = std::lower_bound(mutations.begin(),mutations.end(),-no_categories+2);
	subset.insert(subset.end(),mutations.begin(),upper);

	results.push_back(calculate_frequency_spectrum(subset,n));
	return results;
}

#endif