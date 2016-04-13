/**
 *
 * @author: Marcus Rosti
 * @author: Atallah Hezbor
 */

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

// Array initializer
int **allocate_cells(int n_x, int n_y);

// Given a tour of cities, return an array of all all possible next tours
int ** generate_subproblems(int * tour, int tour_size, int num_of_cities);

void initialize_city_distances(char *filename, int **array, const int num_of_cities);

// Error "handling" -- just kills the process
void die(const char *error) __attribute__ ((const)) __attribute__ ((noreturn));

// Processing workhorses
void master(int **city_dist, const int num_of_cities, const int my_rank, const int nprocs);

void slave(int **city_dist, const int num_of_cities, const int my_rank, const int nprocs);

// Helper function, gets the value of tour
int calculate_tour_distance(int *tour, int tour_size, int **distances) __attribute__((pure));

int main(int argc, char **argv) {
    /////////////////////////////
    // Initialize MPI
    MPI_Init(&argc, &argv);
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (argc < 2)
        die("ERROR: Must provide 2 arguments\n\tmpirun -n {num procs} EXEC"
                    " {number of cities} {distance file}\n");
    const int num_of_cities = atoi(argv[1]);
    char *file_location = argv[2];
    //
    /////////////////////////////

    /////////////////////////////
    // Array Init
    int **cityDistances = allocate_cells(num_of_cities, num_of_cities);
    initialize_city_distances(file_location, cityDistances, num_of_cities);
    //
    /////////////////////////////

    /////////////////////////////
    // Pick a coordination node / or just make it 0
    // TODO refactor this to be distributed
    if (rank == 0) // Make the first processor the master
        master(cityDistances, num_of_cities, rank, nprocs);    
    else // Otherwise their supporting roles
        slave(cityDistances, num_of_cities, rank, nprocs);
    //
    /////////////////////////////

    /////////////////////////////
    // TODO wrap up?
    // Probably need to free things and collect and print my final form

    //
    /////////////////////////////


    MPI_Finalize(); // Close MPI
}

int calculate_tour_distance(int *tour, int tour_size, int **distances) {
    int i;
    int distance = 0;
    // Calculate distance to end of tour
    for (i = 0; i < tour_size - 1; i++) {
        distance += distances[tour[i]][tour[i + 1]];
    }
    // Add distance back to start
    distance += distances[tour[tour_size]][tour[0]];
    return distance;
}

// Allocates and returns a pointer to a 2D array of ints
// TODO could be an error
int **allocate_cells(int num_cols, int num_rows) {
    int **array = (int **) malloc((unsigned int) num_rows * sizeof(int *));
    if (array == NULL) die("Error allocating array!\n");

    array[0] = (int *) malloc((unsigned int) num_rows * (unsigned int) num_cols * sizeof(int));
    if (array[0] == NULL) die("Error allocating array!\n");

    int i;
    for (i = 1; i < num_rows; i++) {    	
        array[i] = array[0] + (i * num_cols);
    }
    
    return array;
}

// Prints the specified error message and then exits
void die(const char *error) {
    printf("%s", error);
    exit(1);
}

void initialize_city_distances(char *filename,              // File pointer to city file
                               int **array,                 // The array to fill
                               const int num_of_cities)     // The number of cities in the file
{
    FILE *fp = fopen(filename, "r");

    int i, j;
    int read_in_number;
    for (i = 0; i < num_of_cities; i++) {
        for (j = 0; j < num_of_cities; j++) {
            fscanf(fp, "%i", &read_in_number);
            array[i][j] = read_in_number;
        }
    }
    fflush(fp);
    fclose(fp);
}

int ** generate_subproblems(int * tour, int tour_size, int num_of_cities) 
{	
	int i, j;
	bool duplicate;
	// there will be at most num_of_cities subproblems
	// each with one more column than the current tour
	// TODO: tighter upper bound?
	int ** subproblems = allocate_cells(num_of_cities, num_of_cities);		
	for (i = 0; i < num_of_cities; i++) {		
		duplicate = false;
		// Check if city already in tour
		for (j = 0; j < tour_size; j++) 
		{			
			if (tour[j] == i) 
			{
				duplicate = true;
				break;
			}
		}
		if (duplicate)		
			continue;
		else 
		{				

			// copy the original tour								
			memcpy((subproblems[i]), tour, tour_size * sizeof(int));			
			// and add the new city to the end
			subproblems[i][tour_size] = i;
		}			
	}
	return subproblems;
}

void master(int **city_dist, const int num_of_cities, const int my_rank, const int nprocs) {
    // Local vars
    int global_lowest_cost = INT32_MAX; // The best path cost
    int *best_path = malloc(num_of_cities * sizeof(int)); // The best path
    best_path[0] = 0;   // Pick a random spot to start
    // So why not 0

    // Generate #nprocs sub problems
    int **subproblems = allocate_cells(nprocs, num_of_cities);

    for (int i = 0; i < nprocs; i++) {

    }



}

/** Algorithm
 * Receive orders from Coordinator
 *
 * Calculate a the cost of the path / find best path given results
 *
 */
void slave(int **city_dist, const int num_of_cities, const int my_rank, const int nprocs) {
    int local_lowest_cost = INT32_MAX;
    int *my_path = malloc(num_of_cities * sizeof(int));

}

