#ifndef RAND_GENERATOR_H
#define RAND_GENERATOR_H

#define LARGE_NUMBER 1000000000
#define LARGE_SEED 161803398
#define MZ 0
#define FAC (1.0 / LARGE_NUMBER);

#include <ctime>
#include <stdlib.h>

#include <vector>
#include <omp.h>

#include <tr1/memory>

class randGenerator;
typedef std::tr1::shared_ptr<randGenerator> spRNG;

using namespace std;

class randGenerator {
public:
	float random();
	int random(int max);
	int random(int min, int max);
	bool mutate(float prob);
	static spRNG getInstance();

private:
	randGenerator(long *dum);


	int _next_;
	int _next_p;
	long ma[56];
	long mj;
	long mk;
	int i;
	int ii;
	int k;
	static vector<spRNG> instances;
};
#endif
