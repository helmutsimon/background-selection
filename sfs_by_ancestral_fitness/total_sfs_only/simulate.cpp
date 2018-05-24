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

void print_frequency_spectrum(Measurement & measurement){
	// std::cout << "Frequency spectrum \n";
	for (auto frequency : measurement.frequency_spectrum)
		std::cout << frequency << ' ';
	std::cout << endl;
}


void print_results(Results & results, bool LAST_ONLY = true){
	if (LAST_ONLY)
		print_frequency_spectrum(results.back());
	else{
		for (auto & measurement:results)
			print_frequency_spectrum(measurement);
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

		//initialize simulation parameters
		double N = strtod(argv[1],NULL);
		double Ud = strtod(argv[2],NULL)/N;
		double Un = strtod(argv[3],NULL)/N;
		double s = -strtod(argv[4],NULL)/N;
		int num_samples = stoi(argv[5],NULL);

		DeltaDFE non_neutral_dfe{s,Ud};
		DeltaDFE neutral_dfe{0,Un};

		Results results;

		results = evolve(random, label_generator, N, non_neutral_dfe, neutral_dfe, num_samples); 
		if (results.size() > 0){
			print_results(results);
		}
		else {
			std::cout<< "ERROR: recorded no frequency spectra.\n";
			return 3;
		}
		return 0;
	}
}
