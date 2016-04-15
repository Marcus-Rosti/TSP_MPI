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
#include <stdint.h>
#include <list.h>


// Array initializer
int **allocate_cells(int n_x, int n_y);

// Given a tour of cities, return an array of all all possible next tours
int fact(const int) __attribute__ ((const));

int ** generate_all_tours_of_depth(const int, const int);

int ** generate_subproblems(int * tour, const int tour_size, const int num_of_cities);

void initialize_city_distances(char *filename, int **array, const int num_of_cities);

// Error "handling" -- just kills the process
void die(const char *error) __attribute__ ((const)) __attribute__ ((noreturn));

// Processing workhorses
void master(int **city_dist, const int num_of_cities, const int my_rank, const int nprocs, const int size_of_work);

void slave(int **city_dist, const int num_of_cities, const int my_rank, const int nprocs, const int size_of_work);

// Helper function, gets the value of tour
int calculate_tour_distance(int *tour, const int tour_size, int **distances, const int num_cities);

void printPath(const int num_of_cities, int *path);

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
        master(cityDistances, num_of_cities, rank, nprocs, 2);
//    else // Otherwise their supporting roles
//        slave(cityDistances, num_of_cities, rank, nprocs, 2);
    //
    /////////////////////////////

    /////////////////////////////
    // TODO wrap up?
    // Probably need to free things and collect and print my final form

    //
    /////////////////////////////


    MPI_Finalize(); // Close MPI
}

int calculate_tour_distance(int *tour, const int tour_size, int **distances, const int num_cities) {
    int i;
    int distance = 0;
    // Calculate distance to end of tour
    for (i = 0; i < tour_size - 1; i++) {
        distance += distances[tour[i]][tour[i + 1]];
    }
    // Add distance back to start
    if(tour_size == num_cities) distance += distances[tour[tour_size-1]][tour[0]];
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

int fact(const int n) {
    if(n<=0) die("Invalid entry");
    if(n==1) return 1;
    return n * fact(n-1);
}

int ** generate_all_tours_of_depth(const int depth, const int num_of_cities) {
    const int num_tours_to_generate = fact(num_of_cities-1)/fact(num_of_cities-depth-1);

    int * tour = malloc((unsigned long) num_of_cities*sizeof(int));

    //////////////////////////////////////////////////////////////////

    int ** depth1 = generate_subproblems(tour,1,num_of_cities);

    printPath(num_of_cities, depth1[2]);

    // TODO this desperately needs refactoring
    //////////////////////////////////////////////////////////////////

    int num_cells = (num_of_cities-1)*(num_of_cities-2);
    int ** depth2 = allocate_cells(num_of_cities,num_cells);

    for(int i = 0; i < num_of_cities-1; i++) {
        int ** iterCities = generate_subproblems(depth1[i],2,num_of_cities);

        for(int j = 0; j < (num_of_cities-2); j++)
            memcpy((depth2[(num_of_cities-2)*i + j]),iterCities[j],(unsigned long)  num_of_cities*sizeof(int));
        free(iterCities);
    }

    //////////////////////////////////////////////////////////////////

    int ** depth3 = allocate_cells(num_of_cities,(num_of_cities-1)*(num_of_cities-2)*(num_of_cities-3));

    for(int i = 0; i < (num_of_cities-1)*(num_of_cities-2); i++) {
        int ** iterCities = generate_subproblems(depth2[i],3,num_of_cities);

        for(int j = 0; j < (num_of_cities-3); j++)
            memcpy((depth3[(num_of_cities-3)*i + j]),iterCities[j],(unsigned long)  num_of_cities*sizeof(int));
        free(iterCities);
    }

    //////////////////////////////////////////////////////////////////

    int ** depth4 = allocate_cells(num_of_cities,(num_of_cities-1)*(num_of_cities-2)*(num_of_cities-3)*(num_of_cities-4));

    for(int i = 0; i < (num_of_cities-1)*(num_of_cities-2)*(num_of_cities-3); i++) {
        int ** iterCities = generate_subproblems(depth3[i],4,num_of_cities);

        for(int j = 0; j < (num_of_cities-4); j++)
            memcpy((depth4[(num_of_cities-4)*i + j]),iterCities[j],(unsigned long)  num_of_cities*sizeof(int));
        free(iterCities);
    }

    //////////////////////////////////////////////////////////////////

    int ** depth5 = allocate_cells(num_of_cities,(num_of_cities-1)*(num_of_cities-2)*(num_of_cities-3)*(num_of_cities-4)*(num_of_cities-5));

    for(int i = 0; i < (num_of_cities-1)*(num_of_cities-2)*(num_of_cities-3)*(num_of_cities-4); i++) {
        int ** iterCities = generate_subproblems(depth4[i],5,num_of_cities);

        for(int j = 0; j < (num_of_cities-5); j++)
            memcpy((depth5[(num_of_cities-5)*i + j]),iterCities[j],(unsigned long) num_of_cities*sizeof(int));
        free(iterCities);
    }

    //////////////////////////////////////////////////////////////////

    return depth5;

}

void printPath(const int num_of_cities, int *path) {
    for(int i = 0; i < num_of_cities; i++)
        printf("%2i  -> ",path[i]);
    printf(" 0\n");
}

/**
 * @arg int * tour         : an array representing the tour we currently generated
 * @arg int tour size      : The size of the current tour
 * @arg int num_of_cities  :
 *
 ****************
 *
 * let num_cities = 5
 * let tour = [0,2,0,0,0]
 * let tour_size = 2
 *
 * returns [[0,2,1,0,0],[0,2,3,0,0],[0,2,4,0,0]]
 *
 ****************
 *
 * let num_cities = 6
 * let tour = [0,1,5,4,0,0]
 * let tour_size = 4
 *
 * returns [[0,1,5,4,2,0],[0,1,5,4,3,0]]
 */
int ** generate_subproblems(int * tour, const int tour_size, const int num_of_cities)
{	
	int i, j;
	bool duplicate;

	// there will be at most num_of_cities subproblems
	// each with one more column than the current tour
    // We only need num_of_cities-tour_size to get the number of new cities that will be added
	int ** subproblems = allocate_cells(num_of_cities, num_of_cities-tour_size);
    int index = 0;
	for (i = 1; i < num_of_cities; i++) { // iterate all of greater than 0 cities
		duplicate = false;
		// Check if city already in tour
		for (j = 1; j < tour_size; j++) // check up to current tour_size to see if the city has been used
		{			
			if (tour[j] == i) 
			{
				duplicate = true;
				break;
			}
		}

		if (duplicate) {
            continue;
        } else {
			// copy the original tour								
			memcpy((subproblems[index]), tour, (unsigned long)  tour_size * sizeof(int));
			// and add the new city to the end
			subproblems[index][tour_size] = i;
            index++;
		}			
	}
	return subproblems;
}

void master(int **city_dist, const int num_of_cities, const int my_rank, const int nprocs, const int size_of_work) {
    // Local vars
    int * first_path = malloc((unsigned long) num_of_cities*sizeof(int));
    for(int i = 1; i < num_of_cities; i++) first_path[i]=i;
    first_path[7] = 6;
    first_path[6] = 7;
    int global_lowest_cost = calculate_tour_distance(first_path,num_of_cities,city_dist,num_of_cities);
    printPath(num_of_cities,first_path);
    printf("Cost %i\n",global_lowest_cost);
    int *best_path = malloc((unsigned long)  num_of_cities * sizeof(int)); // The best path

    int ** work_array = generate_all_tours_of_depth(5,num_of_cities);
    int work_index = 0;


}

/** Algorithmf
 * Receive orders from Coordinator
 *
 * Calculate a the cost of the path / find best path given results
 *
 */
void slave(int **city_dist, const int num_of_cities, const int my_rank, const int nprocs, const int size_of_work) {
    int local_lowest_cost = INT32_MAX;
    int *my_path = malloc((unsigned long) num_of_cities * sizeof(int));


    //MPI_SEND(results);

}



