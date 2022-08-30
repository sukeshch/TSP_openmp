#include "iostream"
#include "randGenerator.h"

vector<spRNG> randGenerator::instances;

using namespace std;

randGenerator::randGenerator(long *dum) {
	// int tid = omp_get_thread_num();
	// cout << "Calling the random generator by thread : " << tid << endl;
    mj = labs( LARGE_SEED - labs(*dum) );
	mj %= LARGE_NUMBER;
	ma[55] = mj;
	mk = 1;

	// Now initialize the rest of the table, in slightly random order,
	// with numbers that are not especially random
	for(i = 1; i <= 54; i++) {
		ii = (21 * i) % 55;
		ma[ii] = mk;
		mk = mj - mk;

		if(mk < MZ)
			mk += LARGE_NUMBER;

		mj = ma[ii];
	}

	// We randomize them by "warming up the generator"
	for(k = 1; k <= 4; k++) {
		for(i = 1; i <= 55; i++) {
			ma[i] -= ma[1 + (i + 30) % 55];

			if(ma[i] < MZ)
				ma[i] += LARGE_NUMBER;
		}
	}

	// Prepare indices for our first generated number
	_next_ = 0;
	_next_p = 31;		// the number 31 is special, don't modify!
	*dum = 1;
}

float randGenerator::random() {
	if(++_next_ == 56)
		_next_ = 1;
	if(++_next_p == 56)
		_next_p = 1;
	mj = ma[_next_] - ma[_next_p];
	if(mj < MZ)
		mj += LARGE_NUMBER;
	ma[_next_] = mj;
	return mj * FAC;
}

bool randGenerator::mutate(float prob) {
	bool flag = false;

	if(random() <= prob)
		flag = true;

	return flag;
}

int randGenerator::random(int min, int max) {
	return min + (int)(random() * (max - min) + 0.5);
}

int randGenerator::random(int max) {
	return this->random(1, max);
}

spRNG randGenerator::getInstance() {
	int id = omp_get_thread_num();

	if( instances.empty() ) {
		static long *tmp = new long[32];
		for(int i = 0; i < 32; i++) {
			int id1 = omp_get_thread_num();
			tmp[i] = (-(unsigned)time(NULL)) * i;
			cout << "instances" << id1 << " " << tmp[i] << " " << &tmp[i] << endl;
			instances.push_back( spRNG(new randGenerator( &tmp[i] )) );
		}

	}

	return instances[id];
}
