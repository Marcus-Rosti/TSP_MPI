/**
 * @author: Marcus Rosti
 * @author: Atallah Hazbor
 */

#include <stdio.h>
#include <mpi.h>

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv); // Initialize MPI

    printf("Hello World\n");

    MPI_Finalize(); // Close MPI
}
