// Assignment 4/5
// Peter Ryder, Brian Kovacik, Matt Holmes

/* INCLUDES */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <limits.h>
#include <pthread.h>
#include <string.h>

#include "clcg4.h"
#include <mpi.h>

#include <hwi/include/bqc/A2_inlines.h>

#define CELL_TYPE unsigned int
#define MULTIPLIER 1000 /* will generate ints between 0 and multiplier */

/* FUNCTION DEFINITIONS */
void generate_matrix();
void print_matrix();
void compute_transpose();
void* sum(void* args);
void write_single_file();
void write_multiple_files();
void cleanup();

/* GLOBAL VARS */
unsigned long long start_cycle_time = 0;
unsigned long long end_cycle_time_compute = 0;
unsigned long long end_cycle_time_output = 0;
unsigned long long total_cycle_time_compute = 0;
unsigned long long total_cycle_time_output = 0;

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

    //Begin counting cycle time
    start_cycle_time = GetTimeBase();
    
	MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &g_commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &g_my_rank);

    g_ranks = g_commsize;
    g_rows_per_rank = g_matrix_size / g_commsize;

#if DEBUG
    if (g_my_rank == 0)
		printf("Running simulation with %d ranks and %d threads per rank\n", g_ranks, g_threads_per_rank);
#endif

    InitDefault();

#if DEBUG
    printf("Rank %d of %d has been started and a first Random Value of %lf\n", 
        g_my_rank, g_commsize, GenVal(g_my_rank));
#endif

    MPI_Barrier(MPI_COMM_WORLD);

#if DEBUG
    if (g_my_rank == 0)
    	printf("Rows per rank: %d\n", g_rows_per_rank);
#endif

	/* generate part of matrix for this rank */
	generate_matrix();

	MPI_Barrier(MPI_COMM_WORLD);

	/* compute transpose */
	compute_transpose();
    
    MPI_Barrier(MPI_COMM_WORLD);

    if (g_my_rank == 0)
        printf("Finished Compute Transpose\n");

    pthread_t* threads;
    int** ends = calloc(g_threads_per_rank, sizeof(int*));
    unsigned int i;
    for (i = 0; i < g_threads_per_rank; i++) {
        ends[i] = calloc(2, sizeof(int));
    }

    threads = calloc(g_rows_per_rank-1, sizeof(pthread_t));
    ends[0][0] = 0; ends[0][1] = g_rows_per_rank/g_threads_per_rank;
    sum(ends[0]);
    for (i = 1; i < g_threads_per_rank; i++) {
        ends[i][0] = g_rows_per_rank/g_threads_per_rank*i; ends[i][1] = g_rows_per_rank/g_threads_per_rank*(i+1);
        pthread_create(threads+i-1, NULL, &sum, ends[i]);
    }
    for (i = 1; i < g_threads_per_rank; i++) {
        pthread_join(threads[i-1], NULL);
    }

	MPI_Barrier(MPI_COMM_WORLD);

    if (g_my_rank == 0)
        end_cycle_time_compute = GetTimeBase();
    
	cleanup();
#if DEBUG
	printf("Matrix %d cleaned up\n", g_my_rank);
#endif

    MPI_Finalize();
    
// stop keeping time, and get the total cycle time
    if (g_my_rank == 0) {
        end_cycle_time_output = GetTimeBase();

        total_cycle_time_output = end_cycle_time_output - end_cycle_time_compute;
        total_cycle_time_compute = end_cycle_time_compute - start_cycle_time;

        printf("Compute Cycle Time: %llu\n", total_cycle_time_compute);
        printf("Output Cycle Time: %llu\n", total_cycle_time_output);        
    }

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
        
        //send row data to other ranks
        for(int rank_iter = 0; rank_iter < g_commsize; rank_iter++)
        {
            MPI_Request request;
            MPI_Isend(&(g_matrix[row][rank_iter*g_rows_per_rank]), g_rows_per_rank, MPI_UNSIGNED, rank_iter, g_my_rank*row+row, MPI_COMM_WORLD, &request);
            MPI_Request_free(&request);
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
    CELL_TYPE* temporary_row = calloc(g_rows_per_rank, sizeof(CELL_TYPE));;
    for(unsigned int row = 0; row < g_rows_per_rank; row++)
    {
        g_matrix_t[row] = calloc(g_matrix_size, sizeof(CELL_TYPE));
    }

    for (int rank_iter = 0; rank_iter < g_commsize; rank_iter++)
    {
        for (CELL_TYPE row = 0; row < g_rows_per_rank; row++) {
            MPI_Request request;
            MPI_Status status;
            MPI_Irecv(temporary_row, g_rows_per_rank, MPI_UNSIGNED, rank_iter, rank_iter*row+row, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status);

            for (CELL_TYPE col = 0; col < g_rows_per_rank; col++) {
                g_matrix_t[col][rank_iter*g_rows_per_rank + row] = temporary_row[col];
            }
        }
    }
    free(temporary_row);
}

void* sum(void* args) {
    int ends[2]; ends[0]=*((int*)args); ends[1]=*((int*)args+1);
#if DEBUG
    fprintf(stderr, "I am thread %lu, responsible for rows %i to %i.\n", pthread_self(), ends[0], ends[1]);
#endif
	for (CELL_TYPE row = ends[0]; row < ends[1]; row++) {
		for (CELL_TYPE col = 0; col < g_matrix_size; col++) {
            g_matrix[row][col] += g_matrix_t[row][col];
		}
	}
}

void cleanup() {
	for (CELL_TYPE row = 0; row < g_rows_per_rank; row++) {
		free(g_matrix[row]);
	}
	free(g_matrix);
}

/*
 *
 *
//1 single file for all MPI ranks, with 8MB block boundaries between rank using MPI File write at all collective IO operation.
void write_single_file()
{
    MPI_Status file_status;
    MPI_File output_file;
    char* filename = "test_out.txt";
    
    unsigned int err = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &output_file);
    if (err != MPI_SUCCESS)
    {
        printf("MPI_FILE_OPEN ERROR PANIC\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int rank_offset = g_my_rank*(g_matrix_size*g_rows_per_rank*sizeof(CELL_TYPE));
    for(int i = 0; i < g_rows_per_rank; i++)
    {
        int offset = i*(g_matrix_size*sizeof(CELL_TYPE));
        err = MPI_File_write_at_all(output_file, rank_offset + offset, g_matrix[i], g_matrix_size, MPI_UNSIGNED, &file_status);
        if (err != MPI_SUCCESS)
        {        
            printf("MPI_FILE_WRITE ERROR PANIC\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    
    err = MPI_File_close(&output_file);
    if (err != MPI_SUCCESS)
    {
        printf("MPI_FILE_CLOSE ERROR PANIC\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}



//# ranks share the same file using MPI File write at (non-collective write call) with 8MB block boundaries between ranks.
void write_multiple_files()
{
    MPI_Status file_status;
    MPI_File output_file;
    MPI_Comm comm_file;
    int file_rank;
    int file_comm_size;
    int split_num = g_my_rank / g_ranks_write_per_file;
    
    MPI_Comm_split(MPI_COMM_WORLD, split_num, g_my_rank, &comm_file);
    MPI_Comm_rank(comm_file, &file_rank);
    MPI_Comm_size(comm_file, &file_comm_size);
    
    char filename[30];
    sprintf(filename, "output%d", split_num);
    
    
    unsigned int err = MPI_File_open(comm_file, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &output_file);
    if (err != MPI_SUCCESS)
    {
        printf("MPI_FILE_OPEN ERROR PANIC\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int rank_offset = g_my_rank*(g_matrix_size*g_rows_per_rank*sizeof(CELL_TYPE));
    for(int i = 0; i < g_rows_per_rank; i++)
    {
        int offset = i*(g_matrix_size*sizeof(CELL_TYPE));
        err = MPI_File_write_at(output_file, rank_offset + offset, g_matrix[i], g_matrix_size, MPI_UNSIGNED, &file_status);
        if (err != MPI_SUCCESS)
        {        
            printf("MPI_FILE_WRITE ERROR PANIC\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    
    err = MPI_File_close(&output_file);
    if (err != MPI_SUCCESS)
    {
        printf("MPI_FILE_CLOSE ERROR PANIC\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
}

*/