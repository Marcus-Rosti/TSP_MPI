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
#include <time.h>

// Array initializer
int **allocate_cells(int n_x, int n_y);

// Given a tour of cities, return an array of all all possible next tours
int fact(const int) __attribute__ ((const)) __attribute__((unused));

int ** generate_all_tours_of_depth(const int, const int);

int ** generate_subproblems(int * tour, const int tour_size, const int num_of_cities);

void initialize_city_distances(char *filename, int **array, const int num_of_cities);

// Error "handling" -- just kills the process
void die(const char *error) __attribute__ ((const)) __attribute__ ((noreturn));

// Processing workhorses
void master(int **city_dist, const int num_of_cities, const int my_rank, const int nprocs, const int size_of_work);

void slave(const int **city_dist, const int num_of_cities, const int my_rank, const int nprocs, const int size_of_work);

// Helper function, gets the value of tour
int calculate_tour_distance(int *tour, const int tour_size, const int **distances, const int num_cities);

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
    const int **cityDistances = allocate_cells(num_of_cities, num_of_cities);
    initialize_city_distances(file_location, cityDistances, num_of_cities);
    //
    /////////////////////////////
    int initial_size = ((num_of_cities-1)*(num_of_cities-2)*(num_of_cities-3)*(num_of_cities-4)*(num_of_cities-5));    
    int work_per_proc = initial_size / nprocs;
    /////////////////////////////
    // Pick a coordination node / or just make it 0
    // TODO refactor this to be distributed
    time_t start_time = time(NULL);
    if (rank == 0) // Make the first processor the master
        master(cityDistances, num_of_cities, rank, nprocs, 10);
    else // Otherwise their supporting roles
        slave(cityDistances, num_of_cities, rank, nprocs, 10);
    printf("Time to calc: %li\n", time(NULL) - start_time);
    /////////////////////////////

    /////////////////////////////
    // TODO wrap up?
    // Probably need to free things and collect and print my final form

    //
    /////////////////////////////


    MPI_Finalize(); // Close MPI
}

int calculate_full_tour_distance(int *tour, const int **distances, const int num_cities) {
    int i;
    int distance = distances[0][tour[1]];
    // Calculate distance to end of tour
    for (i = 1; i < num_cities-1; i++) {
        if(tour[i]   <= 0)          return INT32_MAX;
        if(tour[i]   >= num_cities) return INT32_MAX;
        if(tour[i+1] <= 0)          return INT32_MAX;
        if(tour[i+1] >= num_cities) return INT32_MAX;
        distance += distances[tour[i]][tour[i + 1]];
    }

    // Add distance back to start
    distance += distances[tour[num_cities-1]][tour[0]];
    return distance;
}

int calculate_tour_distance(int *tour, const int tour_size, const int **distances, const int num_cities) {
    if(tour_size == num_cities) {
        return calculate_full_tour_distance(tour, distances, num_cities);
    } else {
        int i;
        int distance = distances[0][tour[1]];
        // Calculate distance to end of tour
        for (i = 1; i < tour_size - 1; i++) {
            if(tour[i] == 0) return INT32_MAX;
            if(tour[i+1] == 0) return INT32_MAX;
            distance += distances[tour[i]][tour[i + 1]];
        }
        // Add distance back to start
        return distance;
    }
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
    //const int num_tours_to_generate = fact(num_of_cities-1)/fact(num_of_cities-depth-1);

    int * tour = malloc((unsigned long) num_of_cities*sizeof(int));
    tour[0] = 0;

    //////////////////////////////////////////////////////////////////

    int ** depth1 = generate_subproblems(tour,1,num_of_cities);

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
    free(depth1);
    free(depth2);
    free(depth3);
    free(depth4);
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

bool valid_path(int * tour, const int num_of_cities) {
    int * seen = calloc(num_of_cities,sizeof(int));
    for(int i = 0; i < num_of_cities; i++) {
        seen[tour[i]] = 1;
    }
    for(int i = 0; i < num_of_cities; i++) {
        if(seen[i] == 0) return false;
    }
    return true;
}

void master(int **city_dist, const int num_of_cities, const int my_rank, const int nprocs, const int size_of_work) {
    // Local vars
    MPI_Status stat;
    MPI_Request req;

    int *first_path = malloc((unsigned long) num_of_cities * sizeof(int));
    int *received_path = malloc((unsigned long) num_of_cities * sizeof(int));
    int received_value;
    for (int i = 0; i < num_of_cities; i++) first_path[i] = i;
    first_path[7] = 6;
    first_path[6] = 7;

    int global_lowest_cost = calculate_tour_distance(first_path, num_of_cities, city_dist, num_of_cities);
    int *best_path = malloc((unsigned long) num_of_cities * sizeof(int)); // The best path

    int **work_array = generate_all_tours_of_depth(5, num_of_cities);
    int work_index = 0;
    int initial_size = ((num_of_cities - 1) * (num_of_cities - 2) * (num_of_cities - 3) * (num_of_cities - 4) *
                        (num_of_cities - 5));
    //int proc_index = 0;
    int kill_signal = 0;


    while (work_index != initial_size) {
        // Send work and bound to each process
        // TODO: handle out of bounds when work not evenly divided
        for (int i = 1; i < nprocs; i++) {
            MPI_Isend(&kill_signal, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &req);
            MPI_Isend(work_array[work_index], num_of_cities * size_of_work, MPI_INT, i, 0, MPI_COMM_WORLD, &req);
            MPI_Isend(&global_lowest_cost, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &req);
            work_index += size_of_work;
        }

        // Receive best path from each process
        for (int j = 1; j < nprocs; j++) {
            MPI_Recv(received_path, num_of_cities, MPI_INT, j, 0, MPI_COMM_WORLD, &stat);
            MPI_Recv(&received_value, 1, MPI_INT, j, 0, MPI_COMM_WORLD, &stat);
            if (received_value < global_lowest_cost) {
                global_lowest_cost = received_value;
                // swap best and received to not lose track of pointers
                memcpy(best_path,received_path,(unsigned long) num_of_cities * sizeof(int));
            }
        }
    }
    // Send kill signal to other processes
    kill_signal = -1;
    for (int i = 1; i < nprocs; i++) {
        MPI_Isend(&kill_signal, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &req);
    }

}

/**
 * @arg int * tour          : partial path to search under
 * @arg int num_of_cities   : total number of cities possible
 * @arg int ** city_dist    : distance matrix
 * @arg int current_size    : current tour length
 * @arg int local_best      : best length according to the caller
 */
int * dfs(int * tour, const int num_of_cities, const int **city_dist, const int current_size, int local_best) {
    int * my_best_path = malloc((unsigned long) num_of_cities * sizeof(int));
    //printf("%i/%i\n",current_size,num_of_cities);
    // Check to see if we've reached the max size, just return the last tour possible
    if(current_size == num_of_cities) {
        //printPath(num_of_cities,tour);
        memcpy(my_best_path, tour, (unsigned long) num_of_cities * sizeof(int));
        //printPath(num_of_cities,my_best_path);
    } else {
        int my_best = local_best;

        memcpy(my_best_path, tour, (unsigned long) num_of_cities * sizeof(int)); // just assume best path is the current tour

        // If we're still alive, generate all of the subproblems

        int **subproblems;
        subproblems = generate_subproblems(my_best_path, current_size, num_of_cities);
        const int num_subproblems = num_of_cities - current_size; // the number of possible subs remaining
        //printf("Generating subprobs: %3i for\t",current_size);
        //printPath(num_of_cities,tour);
        // Now loop over all of the subproblems
        for (int i = 0; i < num_subproblems; i++) {
            // calculate the sub path cost for path at i
            const int sub_path_cost = calculate_tour_distance(subproblems[i], current_size, city_dist, num_of_cities);

            if (sub_path_cost > my_best) continue; // if it exceeds the cost continue

            // otherwise get the best path from my subproblem
            int *path;
            path = dfs(subproblems[i], num_of_cities, city_dist, current_size+1, my_best);

            // if the best of my subproblem is better than local best
            int tempCost = calculate_full_tour_distance(path, city_dist, num_of_cities);
            if (tempCost < my_best && valid_path(path,num_of_cities)) {
                memcpy(my_best_path, path, (unsigned long) num_of_cities * sizeof(int)); // copy it into best path
                my_best = tempCost;
                printf("New best path: %i\t",my_best);
                printPath(num_of_cities,my_best_path);
            }
            free(path);
        }
                    
        free(subproblems[0]);
        free(subproblems);

    }
    //printf("Last dfs path %i --- \t",current_size); printPath(num_of_cities,my_best_path);
    return my_best_path;
}

/** Algorithm
 * Receive orders from Coordinator
 *
 * Calculate a the cost of the path / find best path given results
 *
 */
void slave(const int **city_dist, const int num_of_cities, const int my_rank, const int nprocs, const int size_of_work) {
    MPI_Status stat;
    int local_lowest_cost = INT32_MAX;
    int stay_alive = 1;

    int *local_best_path = malloc((unsigned long) num_of_cities * sizeof(int));
    int ** my_work = allocate_cells(num_of_cities, size_of_work);
    while (true) {
        MPI_Recv(&stay_alive, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
        if (stay_alive == -1) {
            printf("Kill the %ith slave\n",my_rank);
            return;
        }

        // Receive work and bound from master
        MPI_Recv(my_work[0], num_of_cities * size_of_work, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
        MPI_Recv(&local_lowest_cost, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);

        for(int i = 0; i < size_of_work; i++) {
            int * best_dfs_path;
            best_dfs_path = dfs(my_work[i], num_of_cities, city_dist, 6, local_lowest_cost);

            int best_dfs_cost = calculate_full_tour_distance(best_dfs_path, city_dist, num_of_cities);

            if(local_lowest_cost > best_dfs_cost) {
                local_lowest_cost = best_dfs_cost;
                memcpy(local_best_path, best_dfs_path, (unsigned long) num_of_cities * sizeof(int));
                printf("Lowest cost: %i \t",local_lowest_cost);
                printPath(num_of_cities,local_best_path);
            }
            free(best_dfs_path);
        } // Finished all work
        // Send best path and distance to master
        MPI_Send(local_best_path, num_of_cities, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&local_lowest_cost, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }


}
