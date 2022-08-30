**Travelling sales man problem in Assignment 2 of COL380: Introduction to Parallel and distributed programming**


## Running the file

1. g++ main.cpp city.h chromosome.h population.h -fopenmp
2. ./a.out <INPUTFILE> <NUM_THREADS> <POPULATION_SIZE> <NUM_ITERATIONS>

## Design Decisions
1. A random generator file is used to creating random numbers and shared ptr is used over threads to avoid getting same random for diffrent threads.
2. Class for chromosome, city, population are created for modularity.
3. Cities are maintained in an vector and chromosome pattern is maintained by the location of city in this array.
4. Parent for recombination is selected randomly one from elite chromosomes and other from entire population.
5. Both the children produced after pmx recombination are selected.
6. Fittest are selected as top 10 percentage according to distance/fitness.
7. Termination is done when there is no progress in the fittest chromosome for NUM_ITERATIONS.


## Parallelization strategy
1. Code is written such that no two threads write to the same location simultaneously.
2. Parallelization is done using '#pragma omp parallel' for all cases.
3. Parallelized Initial population generation, for every generation - recombinaiton, mutation.

## Load-balancing strategy
1. As every iteration have almost same work load #pragma omp parallel is used for all the parallel parts.
