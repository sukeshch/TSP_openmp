// #include "city.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <fstream>
#include <omp.h>
#include "population.h"
#include "randGenerator.h"
#include <random>

using namespace std;

extern int NUM_THREADS;
extern int DIMENSION;
int POPULATION_SIZE;
int NUM_FITTEST;
int MAX_ITERATIONS;

extern int DIMENSION;
// extern int NUM_FITTEST;
// extern int POPULATION_SIZE;

extern city *cities;
extern double **distances;

extern double eucledian_distance(int id1, int id2);
extern void calculate_distance_between_cities(int dimension);

float PX = 0.8;         // cross over probability
float PMutation = 0.1;  // Mutation probability
population *current_population;
population *new_population;
void selection();
std::pair<chromosome,chromosome> pmx_crossover(chromosome parent1, chromosome parent2);
chromosome mutate(chromosome mychromosome);
double randomInclusive(double max);
double randomExclusive(double max);

int main(int argc, char** argv){
	string INPUT_FNAME = argv[1];
	NUM_THREADS = atoi(argv[2]);
	POPULATION_SIZE = atoi(argv[3]);
	MAX_ITERATIONS = atoi(argv[4]);
	NUM_FITTEST = 10*POPULATION_SIZE/100;

	ifstream fInput;
	fInput.open(INPUT_FNAME);
	fInput >> DIMENSION;
	cities = new city [DIMENSION];
	/* DEBUG */
	cout << INPUT_FNAME << " T" << NUM_THREADS << " D" << DIMENSION << " F"<< NUM_FITTEST << endl;
	cout << "Cities" << endl;
	
	for (int i = 0; i < DIMENSION; ++i)
	{
		int temp_id;
		double temp_x, temp_y;
		fInput >> temp_id >> temp_x >> temp_y;
		cities[i].id = temp_id;
		cities[i].x = temp_x;
		cities[i].y = temp_y;
		
		/* DEBUG */
		// cout << cities[i].id << " " << cities[i].x << " " << cities[i].y << endl;
	}
	fInput.close();
	distances = new double* [DIMENSION];
	for (int i = 0; i < DIMENSION; ++i)
	{
		distances[i] = new double [DIMENSION];
	}
	calculate_distance_between_cities(DIMENSION);
	
	/* DEBUG */
	// for (int i = 0; i < DIMENSION; ++i)
	// {
	// 	for (int j = 0; j < DIMENSION; ++j)
	// 	{
	// 		// cout << i << " " << j << " " << distances[i][j] << endl;
	// 	}
	// }
	cout << "Start of TSP algo \n";
	double tsp_start = omp_get_wtime();
	current_population = new population(POPULATION_SIZE, 1); // INITIALIZED
	new_population = new population(POPULATION_SIZE, 0); // INITIALIZED
	int num_iterations = 0;
	int no_change = 0;
	double old_best_fitness = 0;
	double new_best_fitness = 0;
    while(no_change < MAX_ITERATIONS){
    	cout << "#" << num_iterations << " ";
        selection(); // breed new population
        srand(rand());
        // cout << "Crossover Starts : ";

        /*****************************************************************/
        #pragma omp parallel for num_threads(NUM_THREADS)
        for (int i = NUM_FITTEST; i < POPULATION_SIZE; i+=2)
        {
            std::pair <chromosome, chromosome> pmx_child;
            srand(rand());
            int id1, id2;
            id1 = rand()%POPULATION_SIZE; // Crossover on fittest chromosomes ??
            id2 = rand()%NUM_FITTEST;
            while(id1==id2 || id1 > NUM_FITTEST){
                id1 = rand()%NUM_FITTEST;
            }
            // cout << id1 << " " << id2 << "\n";
            chromosome parent1 = new_population->get_chromosome(id1);
            chromosome parent2 = new_population->get_chromosome(id2);
            // new_population->get_chromosome(0).print_genes();
            if(randomInclusive(1) <= PX){
                // perform cross over
                pmx_child = pmx_crossover(parent1, parent2);
            }
            else{
                // parent is child
                pmx_child = std::make_pair(parent1, parent2);
            }
            current_population->copy_from(mutate(pmx_child.first), i);
            current_population->copy_from(mutate(pmx_child.second), i+1);
        }	
        // cout << "ends \n";
        new_best_fitness = current_population->get_fittest().get_distance();
        if(new_best_fitness==old_best_fitness){
        	no_change++;
        }
        else{
        	no_change=0;
        }
    	old_best_fitness = new_best_fitness;
    	num_iterations++;
    }
	double tsp_end = omp_get_wtime();
    cout << "\n End of GA \n";
    chromosome result = current_population->get_fittest();
    cout << "Result: " << result.get_distance() <<  endl;
    cout << "Total time taken : " << (tsp_end - tsp_start) << endl;

    ofstream fOutput;
    fOutput.open(INPUT_FNAME+"output.txt");
	fOutput << "DIMENSION: " << DIMENSION << endl;
	fOutput << "TOUR LENGTH: " << POPULATION_SIZE << endl;
	fOutput << "TOUR SECTION: \n";
	for (int i = 0; i < DIMENSION; ++i)
	{
		fOutput << result.get_gene(i+1) << endl;
	}
	fOutput << "-1";
	fOutput.close();
	return 0;
}

double randomInclusive(double max)
{
    /* Generate random number r, 0.0 <= r <= max */
    //return ((double)rand() / (double)RAND_MAX * max);
    return ((double)rand() * max) / (double)RAND_MAX;
}

double randomExclusive(double max)
{
    /* Generate random number r, 0.0 <= r < max */
    //return ((double)rand() / ((double)RAND_MAX + 1) * max);
    return ((double)rand() * max) / ((double)RAND_MAX + 1);
}

void selection(){
	// cout << "Inside Selection : \n";
	current_population->sort_population();
    current_population->print_fitness(1);
    #pragma omp parallel for num_threads(NUM_THREADS)
	for (int i = 0; i < NUM_FITTEST; ++i)
	{
		chromosome temp = current_population->get_chromosome(i);
		new_population->copy_from(temp, i);
		// 	cout << "copy Test ";
	    // 	current_population->get_chromosome(i).print_genes();
	    // 	new_population->get_chromosome(i).print_genes();
	}
    // cout << "Printing distances : ";
    // new_population->print_fitness(10);
	// // delete[] current_population;

	population *temp;
	temp = current_population;
	current_population = new_population;
	new_population = temp;
}

std::pair<chromosome,chromosome> pmx_crossover(chromosome parent1, chromosome parent2){ // cross over
	// cout << "Crossover ";
    int startPos = (rand()%DIMENSION) + 1;
    int endPos = (rand()%DIMENSION) + 1;
    while(startPos==endPos){
        endPos = (rand()%DIMENSION) + 1;
    }
    chromosome child1,child2;
    int child1contains[DIMENSION+1],child2contains[DIMENSION+1];
    // #pragma omp parallel for
    for (int i = 1; i <= DIMENSION; ++i)
    {
        child1contains[i] = -1;
        child2contains[i] = -1;
    }
    #pragma omp parallel for num_threads(NUM_THREADS)
    for (int i = 1; i <= DIMENSION; i++) {
        // If our start position is less than the end position
        if (startPos < endPos && i > startPos && i < endPos) {
            int p1 = parent1.get_gene(i);
            int p2 = parent2.get_gene(i);
            child1.set_gene(i, p2);
            child1contains[p2] = p1;

            child2.set_gene(i, p1);
            child2contains[p1] = p2;
        } // If our start position is larger
        else if (startPos > endPos) {
            if (i < startPos && i > endPos) {
                int p1 = parent1.get_gene(i);
	            int p2 = parent2.get_gene(i);
	            child1.set_gene(i, p2);
	            child1contains[p2] = p1;

	            child2.set_gene(i, p1);
	            child2contains[p1] = p2;
            }
        }
    }

    for (int i = 1; i <= DIMENSION; i++) {
        if (startPos < endPos) {
            if (i <= startPos || i >= endPos) {
                if(child1contains[parent1.get_gene(i)] != -1){
                    int temp=parent1.get_gene(i);
                    while(child1contains[temp] != -1){
                        // Unique check of genes in a chromosome
                        temp = child1contains[temp];
                    }
                    child1.set_gene(i, temp);
                    child1contains[child1.get_gene(i)] = parent1.get_gene(i);
                }
                else{ // Not in contains so retain 
                    child1.set_gene(i, parent1.get_gene(i));
                    child1contains[parent1.get_gene(i)] = child1.get_gene(i);
                }

                if(child2contains[parent2.get_gene(i)] != -1){
                    int temp=parent2.get_gene(i);
                    while(child2contains[temp] != -1){
                        // Unique check of genes in a chromosome
                        // cout << temp << child2contains[temp] << endl;
                        temp = child2contains[temp];
                    }
                    child2.set_gene(i, temp);
                    child2contains[child2.get_gene(i)] = parent2.get_gene(i);
                }
                else{ // Not in contains so retain 
                    child2.set_gene(i, parent2.get_gene(i));
                    child2contains[parent2.get_gene(i)] = child2.get_gene(i);
                }
            }

        } // If our start position is larger
        else if (startPos > endPos) {

            if ((i >= startPos || i <= endPos)) {
                
                if(child1contains[parent1.get_gene(i)] != -1){
                    int temp=parent1.get_gene(i);
                    while(child1contains[temp] != -1){
                        // Unique check of genes in a chromosome
                        temp = child1contains[temp];
                    }
                    child1.set_gene(i, temp);
                    child1contains[child1.get_gene(i)] = parent1.get_gene(i);
                }
                else{ // Not in contains so retain 
                    child1.set_gene(i, parent1.get_gene(i));
                    child1contains[parent1.get_gene(i)] = child1.get_gene(i);
                }

                if(child2contains[parent2.get_gene(i)] != -1){
                    int temp=parent2.get_gene(i);
                    while(child2contains[temp] != -1){
                        // Unique check of genes in a chromosome
                        temp = child2contains[temp];
                    }
                    child2.set_gene(i, temp);
                    child2contains[child2.get_gene(i)] = parent2.get_gene(i);
                }
                else{ // Not in contains so retain 
                    child2.set_gene(i, parent2.get_gene(i));
                    child2contains[parent2.get_gene(i)] = child2.get_gene(i);
                }
            }
        }
    }
    // cout << "parent 1 :" ;
    // parent1.print_genes();
    // cout << "parent 2 :" ;
    // parent2.print_genes();
    // cout << "child 1 :" ;
    // child1.print_genes();
    // cout << "child 2 :" ;
    // child2.print_genes();
    return std::make_pair(child1, child2);
}

chromosome mutate(chromosome mychromosome){
	// cout << "Mutate \n";
	// mychromosome.print_genes();
    for(int tourPos1=1; tourPos1 <= DIMENSION; tourPos1++){
        // Apply mutation rate
        if(randomInclusive(1) < PMutation){
            // Get a second random position in the tour
            int tourPos2 = (rand()%DIMENSION)+1;
            while(tourPos2==tourPos1){
            	tourPos2 = (rand()%DIMENSION)+1;
            }

            // Get the cities at target position in tour
            int city1 = mychromosome.get_gene(tourPos1);
            int city2 = mychromosome.get_gene(tourPos2);

            // Swap them around
            mychromosome.set_gene(tourPos2, city1);
            mychromosome.set_gene(tourPos1, city2);
		    // mychromosome.print_genes();
        }
    }
    return mychromosome;
}

