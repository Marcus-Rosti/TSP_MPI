/**
 * @author: Marcus Rosti
 * @author: Atallah Hezbor
 */

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

int calculate_tour_distance(int * tour, int tour_size, int ** distances);
int **allocate_cells(int n_x, int n_y);
void die(const char *error) __attribute__ ((const)) __attribute__ ((noreturn));
void initialize_city_distances(char* filename, int** array, const int num_of_cities);
int ** generate_subproblems(int * tour, int tour_size, int num_of_cities);

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv); // Initialize MPI
    int rank, nthreads;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nthreads);
    if(argc < 2) die("A horrible death\n");
    const int num_of_cities = atoi(argv[1]);
    const char * file_location = argv[2];

    printf("NumCities: %i, File: %s\n",num_of_cities,file_location);

    int** cityDistances = allocate_cells(num_of_cities,num_of_cities);

    initialize_city_distances(file_location,cityDistances,num_of_cities);


    printf("%i\n",cityDistances[10][7]);

    int test_tour[] = {0, 1, 2};        
    printf("An example tour is [%d, %d, %d]\n",test_tour[0], test_tour[1], test_tour[2]);
    int ** new_tours = generate_subproblems(test_tour, 3, 14);
    printf("A new tour is [%d, %d, %d, %d] \n", new_tours[7][0], new_tours[7][1], new_tours[7][2], new_tours[7][3]);

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

// Allocates and returns a pointer to a 2D array of ints
// TODO could be an error
int **allocate_cells(int num_cols, int num_rows)
{
    int **array = (int **) malloc(num_rows * sizeof(int *));
    if (array == NULL) die("Error allocating array!\n");

    array[0] = (int *) malloc(num_rows * num_cols * sizeof(int));
    if (array[0] == NULL) die("Error allocating array!\n");

    int i;
    for (i = 1; i < num_rows; i++) {    	
        array[i] = array[0] + (i * num_cols);
    }
    
    return array;
}


// Prints the specified error message and then exits
void die(const char *error)
{
    printf("%s", error);
    exit(1);
}

void initialize_city_distances(char *filename,              // File pointer to city file
                               int **array,                 // The array to fill
                               const int num_of_cities)     // The number of cities in the file
{
    FILE *fp = fopen(filename,"r");

    int i,j;
    int read_in_number;
    for(i=0; i<num_of_cities; i++)
    {
        for(j = 0; j < num_of_cities; j++)
        {
            fscanf(fp,"%i",&read_in_number);
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
