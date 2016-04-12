/**
 * @author: Marcus Rosti
 * @author: Atallah Hazbor
 */

#include <stdio.h>
#include <mpi.h>





int calculate_tour_distance(int * tour, int tour_size, int ** distances);

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv); // Initialize MPI

    printf("Hello World\n");

    MPI_Finalize(); // Close MPI
}

int calculate_tour_distance(int * tour, int tour_size, int ** distances) {	
	int i;
	int distance = 0;
	// Calculate distance to end of tour
	for (i=0; i < tour_size-1; i++) 
	{
		distance += distances[tour[i]][tour[i+1]];
	}
	// Add distance back to start
	distance += distances[tour[tour_size]][tour[0]];
	return distance;
}
