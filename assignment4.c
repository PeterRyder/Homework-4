// Assignment 4/5
// Peter, Brian, Matt

/* INCLUDES */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <limits.h>

#include <clcg4.h>
#include <mpi.h>

#include <hwi/include/bqc/A2_inlines.h>

#define CELL_TYPE unsigned int
#define MULTIPLIER 1000 /* will generate ints between 0 and multiplier */

/* FUNCTION DEFINITIONS */
void generate_matrix();
void print_matrix();
void compute_transpose();
void write_single_file();
void write_multiple_files();
void cleanup();

/* GLOBAL VARS */
unsigned long long start_cycle_time = 0;
unsigned long long end_cycle_time = 0;
unsigned long long total_cycle_time = 0;

unsigned int g_ranks = -1;
unsigned int g_threads_per_rank = 0;
unsigned int g_matrix_size = 0;
unsigned int g_ranks_write_per_file = 1;

int g_my_rank = -1;
int g_commsize = -1;

unsigned int g_rows_per_rank = 0;

CELL_TYPE **g_matrix = NULL;
CELL_TYPE **g_matrix_t = NULL;

int main(int argc, char* argv[]) {

	if (argc != 3) {
		printf("Wrong arguments\n");
		printf("usage: assignment4 [threads_per_rank] [matrix_size]");
		return -1;
	}

#if DEBUG
	printf("Got %d argument(s)\n", argc);

	for (int i = 0; i < argc; i++) {
		printf("%s\n", argv[i]);
	}
#endif

	g_threads_per_rank = atoi(argv[1]);
	g_matrix_size = atoi(argv[2]);

	if (g_matrix_size <= 0 || g_ranks <= 0) {
		printf("Bad arguments\n");
		return(-1);
	}

	/* create ranks here */
	printf("Initializing MPI\n");
    
    //Begin counting cycle time
    start_cycle_time = GetTimeBase();
    
	MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &g_commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &g_my_rank);

    g_ranks = g_commsize;

#if DEBUG
    if (g_my_rank == 0)
		printf("Running simulation with %d ranks and %d threads per rank\n", g_ranks, g_threads_per_rank);
#endif

    InitDefault();

    printf("Rank %d of %d has been started and a first Random Value of %lf\n", 
	   g_my_rank, g_commsize, GenVal(g_my_rank));

    MPI_Barrier(MPI_COMM_WORLD);

    g_rows_per_rank = g_matrix_size / g_commsize;

#if DEBUG
    if (g_my_rank == 0)
    	printf("Rows per rank: %d\n", g_rows_per_rank);
#endif

	/* generate part of matrix for this rank */
	generate_matrix();

	MPI_Barrier(MPI_COMM_WORLD);

	print_matrix();

	/* compute transpose */
	compute_transpose();

	MPI_Barrier(MPI_COMM_WORLD);

    //PRINT OUTPUT COMMENTED OUT UNTIL I GET A CHANCE TO TEST IT
/*    if( g_commsize <= g_ranks_write_per_file )
    {
        write_single_file(g_my_rank);
    }
    else
    {
        write_multiple_files(g_my_rank);
    }*/
    
	cleanup();
#if DEBUG
	printf("Matrix %d cleaned up\n", g_my_rank);
#endif

    MPI_Finalize();
    
// stop keeping time, and get the total cycle time
    end_cycle_time = GetTimeBase();
    total_cycle_time = end_cycle_time - start_cycle_time;

#if DEBUG
    if (g_my_rank == 0)
    	printf("MPI_Finalize complete\n");
#endif
	
	return 0;
}

void generate_matrix() {
	/* initialize slice of matrix */
	g_matrix = calloc(g_rows_per_rank, sizeof(CELL_TYPE *));

	for (CELL_TYPE row = 0; row < g_rows_per_rank; row++) {
		g_matrix[row] = calloc(g_matrix_size, sizeof(CELL_TYPE));

		for (CELL_TYPE col = 0; col < g_matrix_size; col++) {
			g_matrix[row][col] = GenVal(g_my_rank) * MULTIPLIER;
		}
	}
}

void print_matrix() {
	for (CELL_TYPE row = 0; row < g_rows_per_rank; row++) {
		for (CELL_TYPE col = 0; col < g_matrix_size; col++) {
			printf("%u|", g_matrix[row][col]);
		}
		printf("\n");
	}
}

void compute_transpose() {
	g_matrix_t = calloc(g_rows_per_rank, sizeof(CELL_TYPE *));

	for (CELL_TYPE row = 0; row < g_rows_per_rank; row++) {
		g_matrix_t[row] = calloc(g_matrix_size, sizeof(CELL_TYPE));

		for (CELL_TYPE col = 0; col < g_matrix_size; col++) {
			/* get transposed value */
		}
	}
}

void cleanup() {
	for (CELL_TYPE row = 0; row < g_rows_per_rank; row++) {
		free(g_matrix[row]);
	}
	free(g_matrix);
}

//1 single file for all MPI ranks, with 8MB block boundaries between rank using MPI File write at all collective IO operation.
void write_single_file()
{
    MPI_Status file_status;
    MPI_file output_file;
    char* filename = "test_out.txt" //REPLACE ME
    
    unsigned int err = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_WRONGLY, MPI_INFO_NULL, &output_file);
    if (err != MPI_SUCCESS)
    {
        print_MPI_error(err, "MPI_File_open");
    }
    for(int i = 0; i < g_rows_per_rank; i++)
    {
        err = MPI_File_write_at_all(output_file, g_my_rank, matrix[i], g_matrix_size, MPI_FLOAT, &file_status);
        if (err != MPI_SUCCESS)
        {
            print_MPI_error(err, "MPI_File_write_at_all");
        }
    }
    
    err = MPI_File_close(&output_file);
    if (err != MPI_SUCCESS)
    {
        print_MPI_error(err, "MPI_File_close");
    }
}

//# ranks share the same file using MPI File write at (non-collective write call) with 8MB block boundaries between ranks.
void write_multiple_files()
{
    MPI_Status file_status;
    MPI_file output_file;
    MPI_Comm mpi_comm_file;
    
    //SHARE FLIE
}