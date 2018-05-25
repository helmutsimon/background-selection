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

void print_array(std::valarray<double> & v, int total_repeats){
	for (auto item:v)
		std::cout << item/float(total_repeats) << ' ';
	std::cout << endl;
}


void print_results(Results & results, int total_repeats){
	print_array(results.total_frequency_spectrum,total_repeats);
	for (auto cond_sfs : results.conditional_sfss){
		print_array(cond_sfs,total_repeats);
	}
}

int simulate_single_effect_DFE(int argc, char * argv[]);

int main(int argc, char* argv[])
{	
	return simulate_single_effect_DFE(argc,argv);
}

int simulate_single_effect_DFE(int argc, char* argv[]){
	if (argc < 6)
	{
        std::cout << "usage: " << argv[0] << " N NUd NUn Ns num_samples" << std::endl;
        return 1;
    }
    else{

    	// initialize random number generator and mutation label generator
		Random random = create_random();
		LabelGenerator label_generator{};

		//initialize parameters
		double N = strtod(argv[1],NULL);
		double Ud = strtod(argv[2],NULL)/N;
		double Un = strtod(argv[3],NULL)/N;
		double s = -strtod(argv[4],NULL)/N;
		int num_samples = stoi(argv[5],NULL);

		DeltaDFE non_neutral_dfe{s,Ud};
		DeltaDFE neutral_dfe{0,Un};

		Results results;

		results = evolve(random, label_generator, N, non_neutral_dfe, neutral_dfe, num_samples); 
		print_results(results,num_samples);

		return 0;
	}
}




