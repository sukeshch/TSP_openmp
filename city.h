#include <cmath>

class city{
public:
	int id;
	double x;
	double y;
};

city *cities;
double **distances;
int DIMENSION;
int NUM_THREADS;


double eucledian_distance(int id1, int id2)
{
	double sqx = cities[id1].x - cities[id2].x;
	double sqy = cities[id1].y - cities[id2].y;
	sqx = sqx*sqx;
	sqy = sqy*sqy;
	return sqrt(sqx+sqy);
}

void calculate_distance_between_cities(int DIMENSION)
{
	// #pragma omp parallel for
	for (int i = 0; i < DIMENSION; ++i)
	{
		for (int j = 0; j < i; ++j)
		{
			if(i!=j){
				distances[i][j] = eucledian_distance(i,j);
				distances[j][i] = distances[i][j];
			}
		}
	}
}
