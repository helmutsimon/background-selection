#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>  
using namespace std;

#include "individual.hpp"
#include "dfe.hpp"
#include "evolve.hpp"


void print_mutation(Mutation & mutation){
	std::cout << mutation.label << mutation.fitness << endl;
}


template<class Type>
void print_vector(std::vector<Type> v){
	for (auto & item:v)
		std::cout << item << ' ';
	std::cout<< endl;
}

void print_results(Results & results, int total_repeats){
	std::set<std::pair<int,int>> unique_counts(results.counts.begin(),results.counts.end());
	std::multiset<std::pair<int,int>> all_counts(results.counts.begin(),results.counts.end());
	
	for (auto &item: unique_counts){
		std::cout << item.first << " ";
	}
	std::cout << std::endl;
	for (auto &item: unique_counts){
		std::cout << item.second << " ";
	}
	std::cout << std::endl;
	for (auto &item: unique_counts){
		std::cout << all_counts.count(item) << " ";
	}
	std::cout << std::endl;
}

int simulate_single_effect_DFE(int argc, char * argv[]);

int main(int argc, char* argv[])
{	
	return simulate_single_effect_DFE(argc,argv);
}

int simulate_single_effect_DFE(int argc, char* argv[]){
	if (argc < 7)
	{
        std::cout << "usage: " << argv[0] << " N NUd NUn Ns time_interval num_samples" << std::endl;
        return 1;
    }
    else{

    	// initialize random number generator and mutation labelgenerator
		Random random = create_random();
		LabelGenerator label_generator{};

		//initialize JDFE and population size
		double N = strtod(argv[1],NULL);
		double Ud = strtod(argv[2],NULL)/N;
		double Un = strtod(argv[3],NULL)/N;
		double s = -strtod(argv[4],NULL)/N;
		int time_interval = stoi(argv[5],NULL);
		int num_samples = stoi(argv[6],NULL);


		DeltaDFE non_neutral_dfe{s,Ud};
		DeltaDFE neutral_dfe{0,Un};


		Results results;

		results = evolve(random, label_generator, N, non_neutral_dfe, neutral_dfe, num_samples, time_interval); 
		print_results(results,num_samples);

		return 0;
	}
}




