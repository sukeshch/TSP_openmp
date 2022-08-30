#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <omp.h>
#include "chromosome.h"

int population_size;
extern int DIMENSION;
extern int NUM_THREADS;

bool compare_chromosome(chromosome c1, chromosome c2){
	return c1.get_distance() < c2.get_distance();
}

class population
{
	std::vector<chromosome> chromosomes;
	int current_size;
public:
	population(int size, bool initialize){
		// Initialization of random chromosomes
		population_size = size;
		/* DEBUG */
		// std::cout << "population_size : " << population_size << "\n";
		chromosomes.resize(population_size);
		
		if(initialize){
			#pragma omp parallel for num_threads(NUM_THREADS)
			for (int i = 0; i < population_size; ++i)
			{
				chromosomes[i].generate_chromosome();
			}
		}
		else{
			#pragma omp parallel for num_threads(NUM_THREADS)
			for (int i = 0; i < population_size; ++i)
			{
				chromosomes[i].gene_initialize();
			}	
		}
	}
	~population(){
		// delete[] chromosomes; 
	}

	chromosome get_chromosome(int index){
		return chromosomes.at(index);
	}
	
    void copy_from(population *parent, int index){
    	chromosomes[index].copy_genes_from(parent->get_chromosome(index));
    }
	void copy_from(chromosome parent, int index){
    	chromosomes[index].copy_genes_from(parent);
    }
	void sort_population(){
		// std::cout << "Sort population \n";
		// returns sorted - descending order of chromosomes
		std::sort(chromosomes.begin(), chromosomes.end(), compare_chromosome); 
	}
	
	chromosome get_fittest() {
        chromosome fittest = chromosomes[0];
        // Loop through individuals to find fittest
        for (int i = 1; i < population_size; i++) {
            if (fittest.get_fitness() <= chromosomes[i].get_fitness()) {
                fittest = chromosomes[i];
            }
        }
        return fittest;
    }
    void print_fitness(int count) 
    {
    	for (int i = 0; i < count; ++i)
    	{
    		std::cout << chromosomes[i].get_distance() << " "; 
    	}
    	std::cout << " \n";
    }
};
