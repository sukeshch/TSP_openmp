#include <iostream>
#include <cstdlib>
#include <time.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include "city.h"
extern city *cities;
extern double **distances;
extern int DIMENSION;
extern int NUM_THREADS;

class chromosome
{
	std::vector<int> genes;
	int chromosome_distance = 0;
	double fitness = 0;
	std::vector<int> sequence_vector;
public:
	// chromosome(){
		
	// };
	chromosome(){
		// chromosome_length = length;
		genes.resize(DIMENSION+1);
		for (int i = 0; i < DIMENSION; ++i)
		{
			sequence_vector.push_back(i+1);
		}
	}
	~chromosome(){
		// delete[] genes;
	}


	void generate_chromosome(){
		// Randomly generate chromosomes - Populating the initial set of chromosomes
		srand(rand());
		#pragma omp critical
		{
		random_shuffle(sequence_vector.begin(), sequence_vector.end());
		}
		
		for (int i = 1; i <= DIMENSION; ++i)
		{
			genes[i] = sequence_vector.at(i-1);
		}
	}

	int get_gene(int index){
		if( index <= DIMENSION && index >0 ){
			return genes[index];
		}
		std::cout << "Index out of bounds:"<< index << " get_gene\n";
		return -1; // invalid 
	}	
	void set_gene(int index, int value){
		if( index <= DIMENSION && index>0){
			if(value <= DIMENSION && value > 0 ){
				genes[index] = value;
				return;
			}
			else{
				std::cout << "Value out of bounds:"<< value << " set_gene \n";
			}
		}
		std::cout << "Index out of bounds: "<< index << "get_gene\n";
	}
	double get_fitness() {
        if (fitness == 0) {
            fitness = 1/((double)get_distance());
        }
        return fitness;
    }
    void gene_initialize(){
    	for (int i = 0; i < DIMENSION; ++i)
    	{
    		genes[i+1] = -1;
    	}
    }
    // Gets the total distance of the tour
    int get_distance(){
        if (chromosome_distance == 0) {
            int distance = 0;
            // Loop through our tour's cities
            for (int cityIndex=1; cityIndex <= DIMENSION; cityIndex++) {
                // Get city we're travelling from
                int fromCity = genes[cityIndex];
                // City we're travelling to
                int destinationCity;
                // Check we're not on our tour's last city, if we are set our
                // tour's final destination city to our starting city
                if(cityIndex != DIMENSION){
                    destinationCity = genes[cityIndex+1];
                }
                else{
                    destinationCity = genes[1];
                }
                // Get the distance between the two cities
                distance += distances[fromCity-1][destinationCity-1];
            }
            chromosome_distance = distance;
        }
        return chromosome_distance;
    }
    void copy_genes_from(chromosome parent){
    	// std::cout << "copy_genes_from ";
    	// parent.print_genes();
    	for (int i = 1; i <= DIMENSION; ++i)
    	{
    		int temp = parent.get_gene(i);
    		genes.at(i) = temp;
    		// std::cout << genes.at(i) << " " << parent.get_gene(i) << "\n";
    	}
    	fitness = 0;
    	chromosome_distance=0;
    }

    void print_genes(){
    	for (int i = 1; i <= DIMENSION; ++i)
    	{
    		std::cout << genes[i] << " " ;
    	}
    	std::cout << "\n";
    }
};