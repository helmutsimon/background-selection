#include <vector>
#include <iostream>
#include <cmath>
using namespace std;

#include "stats.hpp"

typedef std::vector<double> FrequencySpectrum;
class FrequencySpectra{
public:
	FrequencySpectrum synonymous;
	FrequencySpectrum nonsynonymous;
};

FrequencySpectra calculate_frequency_spectra(double NUd, double Ns, int n, int trials); 

class FitnessClass{
public:
	double coalescence_rate;
	double mutation_rate;
	int k;
	std::vector<int> lineages;
	FitnessClass(double coalescence_rate=0, double mutation_rate=0, int size_hint = 1): coalescence_rate(coalescence_rate), mutation_rate(mutation_rate) { lineages.reserve(size_hint); }
};

typedef std::vector<FitnessClass> FitnessDistribution;

FitnessDistribution create_fitness_distribution(double Ns, double lambda, int n);
void extend_fitness_distribution(FitnessDistribution & fitness_distribution, int kmax);
void print_fitness_distribution(FitnessDistribution & fitness_distribution);


class Event{
public:
	enum EventType{MUTATION, COALELSCENCE};
	double weight;
	EventType type;
	decltype(FitnessDistribution().begin()) k_ptr; // the fitness class where this event occurs
    bool operator<(const double p) const { return weight < p; } // used for intrusive search		
};

int main(int argc, char* argv[])
{	
	if (argc < 5)
	{
        std::cout << "usage: " << argv[0] << " NUd Ns n num_repeats" << std::endl;
        return 1;
    }
    else{

		// read in parameters
		double NUd = strtod(argv[1],NULL);
		double Ns = strtod(argv[2],NULL);
		int n = stoi(argv[3],NULL);
		int trials = stoi(argv[4],NULL);
		// std::cout << "started" << std::endl;
		// time is measured in units of N generations
		// neutral frequency spectrum should be scaled by NUn
		FrequencySpectra frequency_spectra = calculate_frequency_spectra(NUd, Ns, n, trials);

		//std::cout << "PRINTING FREQUENCY SPECTRUM... \n";
		// for (auto frequency_ptr = frequency_spectra.synonymous.begin()+1; frequency_ptr != frequency_spectra.synonymous.end(); ++frequency_ptr){
		// 	std::cout << *frequency_ptr << ' ';
		// }
		// std::cout << "\n";	

		for (auto frequency_ptr = frequency_spectra.synonymous.begin()+1; frequency_ptr != frequency_spectra.synonymous.end(); ++frequency_ptr){
			std::cout << *frequency_ptr << ' ';
		}
		std::cout << "\n";	
	}

	return 0;
}


FrequencySpectra calculate_frequency_spectra(double NUd, double Ns, int n, int trials){

	// initialize random number generator
	Random random = create_random();

	// initialize fitness distribution
	double lambda = NUd/Ns;

	FrequencySpectra frequency_spectra{FrequencySpectrum(n,0.0),FrequencySpectrum(n,0.0)};
	FitnessDistribution fitness_distribution = create_fitness_distribution(Ns, lambda, n);
	int kmax = fitness_distribution.size()-1;

	for (int current_trial  = 0; current_trial < trials; ++current_trial){
		// randomly allocate n individuals to site frequency spectrum
		// kmax keeps track of the size of the fitness distribution
		for (int i = 0; i < n; ++i){
			int k = static_cast<int>(sample_poisson(random,lambda));
			if (k>kmax){
				extend_fitness_distribution(fitness_distribution, k);	kmax = k;
			}
			fitness_distribution[k].lineages.push_back(1);
		}
		
		auto kmin_ptr = fitness_distribution.begin();
		while(kmin_ptr->lineages.size() == 0){
			++kmin_ptr;
		}

		auto kmax_ptr = fitness_distribution.begin()+kmax;
		while(kmax_ptr->lineages.size() == 0){
			--kmax_ptr;
		}

		std::vector<Event> event_list;

		//print_fitness_distribution(fitness_distribution);

		while(true){
			// compile a list of all possible events
			double total_rate = 0.0;
			event_list.clear();

			for (auto k_ptr = kmin_ptr; k_ptr <= kmax_ptr; ++k_ptr){
				int no_lineages = k_ptr->lineages.size();
				if ((no_lineages>0))
					event_list.push_back(Event{total_rate+=k_ptr->mutation_rate * no_lineages, Event::MUTATION, k_ptr});
				if (no_lineages>1)
					event_list.push_back(Event{total_rate+=k_ptr->coalescence_rate * (no_lineages*(no_lineages-1)/2.0) , Event::COALELSCENCE, k_ptr});
			}

			//pick an event at random
			auto event_ptr = std::lower_bound(event_list.begin(), event_list.end(), total_rate*sample_uniform(random)); 

			double event_time = 1/total_rate;

			//update neutral SFS
			for (auto k_ptr = kmin_ptr; k_ptr <= kmax_ptr; ++k_ptr){
				for(auto i = k_ptr->lineages.begin(); i!=k_ptr->lineages.end(); ++i){
					frequency_spectra.synonymous[*i] += event_time;
				}
			}

			// perform event
			if (event_ptr->type == Event::MUTATION){
				auto k_ptr = event_ptr->k_ptr;
				int mutated_lineage_index = static_cast<int>(sample_uniform(random)*k_ptr->lineages.size());
				//std::cout<< "Mutated lineage "<< mutated_lineage_index << " in class " << k_ptr - fitness_distribution.begin() << "\n";
				frequency_spectra.nonsynonymous[k_ptr->lineages[mutated_lineage_index]]++;
				(k_ptr-1)->lineages.push_back(k_ptr->lineages[mutated_lineage_index]);
				k_ptr->lineages.erase(k_ptr->lineages.begin() + mutated_lineage_index);

				// check whether we need to update lowest and highest populated classes
				if (k_ptr == kmin_ptr){
					--kmin_ptr;
				}

				if (kmin_ptr < fitness_distribution.begin())
				{
					std::cout << "EXCEEDED FITNESS DISTRIBUTION!!!!\n";
				}
				if(kmax_ptr->lineages.size() == 0){
					--kmax_ptr;
				}
			}
			else{
				auto k_ptr = event_ptr->k_ptr;
				int lineage_1_index = static_cast<int>(sample_uniform(random)*k_ptr->lineages.size());
				int lineage_2_index = static_cast<int>(sample_uniform(random)*(k_ptr->lineages.size()-1));
				if (lineage_2_index >= lineage_1_index)
					++lineage_2_index;

				k_ptr->lineages[lineage_1_index] += k_ptr->lineages[lineage_2_index];
				k_ptr->lineages.erase(k_ptr->lineages.begin()+lineage_2_index);
				//std::cout<< "Coalesced lineages "<< lineage_1_index <<" and " << lineage_2_index << " in class " << k_ptr - fitness_distribution.begin() << "\n";

				if (k_ptr->lineages[lineage_1_index] == n){
					k_ptr->lineages.clear();
					break;
				}
			}
		}
	}

	for (int i = 0; i < n; i++){
		frequency_spectra.synonymous[i] = frequency_spectra.synonymous[i]/trials;
		frequency_spectra.nonsynonymous[i] = frequency_spectra.nonsynonymous[i]/trials;
	}
	return frequency_spectra;

}

FitnessDistribution create_fitness_distribution(double Ns, double lambda, int n){
	FitnessDistribution fitness_distribution;
	fitness_distribution.push_back(FitnessClass(exp(lambda),0.0,n));
	fitness_distribution.push_back(FitnessClass(exp(lambda)/lambda, Ns, n));
	extend_fitness_distribution(fitness_distribution,static_cast<int>(1 + 2 * lambda));
	return fitness_distribution;
}

void extend_fitness_distribution(FitnessDistribution & fitness_distribution, int kmax){
	double lambda = fitness_distribution[0].coalescence_rate / fitness_distribution[1].coalescence_rate;
	for (int k = fitness_distribution.size(); k<=kmax; ++k){
		fitness_distribution.push_back(FitnessClass(fitness_distribution.back().coalescence_rate * k / lambda, fitness_distribution.back().mutation_rate * k /(k-1), fitness_distribution.back().lineages.capacity()) );
	}
}

void print_fitness_distribution(FitnessDistribution & fitness_distribution)
{
	for (int i = 0; i < fitness_distribution.size(); ++i){
		std::cout << "Fitness class " << i << ":\t";
		if (fitness_distribution[i].lineages.size() == 0)
			std::cout<< "EMPTY \n"; 
		else{
				for (auto j = fitness_distribution[i].lineages.begin(); j != fitness_distribution[i].lineages.end(); ++j){
					std::cout << *j << ' ';
				}
				std::cout << std::endl;
		}
	}
}