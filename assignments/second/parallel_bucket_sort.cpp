//
// Created by guowenmi on 2019-10-01.
/**
 * parallel bucket sort
 *
 * how to do?
 * it will use master/slave mode.
 *
 * 1. First, Master generates random numbers needed to be sorted. (only in Master)
 * 2. Second, Master scatter numbers to each processors evenly.
 * 3. Each processor divides its own numbers into proper small buckets.
 * 4. Each processor scatters its own small buckets to all other processors
 *    and gather its own proper numbers from all other processors. (by call alltoall() and alltoallv() functions)
 * 5. Each processor sorts its own numbers
 * 6. Master gather sorted numbers from all other processes (large bucket).
 * 7. Master combines these large bucket in proper order.
 *
 * the steps, 3, 4 and 5, will be paralleled.
 *
 */
#include "mpi.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <time.h>
#include <cstdlib>
#include <tr1/random>

#define INF (-1)
using namespace std;

unsigned long number_size;//length of the unsorted data
int curr_rank;//rank of the current process
int processes_number;//the number of processes
int MASTER_RANK = 0; // the master's rank

//generate the unsorted data, range (0, max)
unsigned long *generate_array_with_random(unsigned long max, unsigned long size){
//    std::tr1::default_random_engine e;
//    std::tr1::uniform_int_distribution<unsigned> u(min, max);

    unsigned long *array = new unsigned long [size];
    for (unsigned long i = 0; i < size; i ++) {
    //    array[i] = u(e);
        array[i] = rand() % max + 1;
        cout << "array_random_number = " << array [i] << endl;
    }
    return array;
}

//IncOrder for qsort
int IncOrder(const void *e1, const void *e2)
{
    return (*((unsigned long *)e1)-*((unsigned long *)e2));
}

void display(unsigned long *array, unsigned long size) {
    for(unsigned long i = 0; i<size; i++)
        cout << array[i] << " ";
    cout << endl;
}

void display_int_array(int *array, int size) {
    for(int i = 0; i<size; i++)
        cout << array[i] << " ";
    cout << endl;
}

int main(int argc, char **argv)
{
    if(argc!=2)
    {
        cout<<"Please check your input"<<endl;
        exit(0);
    }
    unsigned long curr_proc_data_size;//size of data on the current process
    unsigned long *curr_proc_data;//data on the current process
    unsigned long *original_data;//the unsorted data
//    long *final_sorted_data;//the sorted data
    unsigned long *pivot_list;//the pivot list
    unsigned long *small_buckets;//the buckets in the current process
    unsigned long *recv_bucket_alltoallv;//the bucket after alltoallv function
    double cost_time;

    MPI_Init(&argc, &argv); //initial

    //Blocks until all processes in the communicator have reached this routine
//    MPI_Barrier(MPI_COMM_WORLD);
    cost_time = - MPI_Wtime();
    MPI_Comm_rank(MPI_COMM_WORLD, &curr_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processes_number);
    number_size = atoll(argv[1]); // the size of numbers need to be sorted, an input parameter
    curr_proc_data_size = number_size / processes_number; // the number of numbers is total_number / processers

    cout << "number_size = " << number_size << ", curr_proc_data_size = " << curr_proc_data_size << endl;


    // step 1, initial the number
    if(curr_rank == MASTER_RANK)
    {
        original_data = generate_array_with_random (number_size, number_size);
    }

    /*
     * these two lines does need if calling scatter.
    MPI_Bcast(&number_size, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&curr_proc_data_size, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    */

    // step 2, scatter evenly to all processes.
    curr_proc_data = new unsigned long[curr_proc_data_size];
//    final_sorted_data=new long[number_size];
    pivot_list=new unsigned long[processes_number];
    MPI_Scatter(original_data, curr_proc_data_size, MPI_LONG, curr_proc_data, curr_proc_data_size, MPI_LONG, MASTER_RANK, MPI_COMM_WORLD);

    for (unsigned long i = 0; i < curr_proc_data_size; i ++){
        cout << "after scatter rank " << curr_rank << " : curr_proc_data [i] = " << curr_proc_data [i] << cost_time<<endl;
    }

    // step 3, each processor loops each number to determine which small bucket they should in.
    //initialize small_buckets
    small_buckets = new unsigned long[number_size];
    for(unsigned long i = 0; i < number_size; i++)
    {
        small_buckets[i] = INF;
    }

    // update small_buckets
    // curr_proc_data --> small_buckets
    int buckets_number = processes_number;
    unsigned long *bucket = (unsigned long*)calloc(buckets_number * curr_proc_data_size, sizeof(long));
    // new unsigned long [buckets_number * curr_proc_data_size]

    //initialize number of items, used to storte the size of numbers in small buckets
    int *nitems = (int*)calloc(buckets_number, sizeof(int));
    unsigned long step = number_size/processes_number;

    for (int i = 0; i < curr_proc_data_size; i++)
    {
        int bktno = floor(curr_proc_data[i]/step);// in which small bucket
        int idx = bktno * curr_proc_data_size + nitems[bktno];// index in the bucket
        bucket[idx] = curr_proc_data[i];
        ++nitems[bktno];
    }

    // step 4, each processor scatter its numbers to proper processors and gather its own proper numbers from others
    // firstly, need to let all processores know how many numbers should recv from each processor
    int* recv_count_alltoallv = (int*)calloc(buckets_number, sizeof(int));
 //   int send_count, recv_count = 1;
    MPI_Alltoall(nitems, 1, MPI_INT, recv_count_alltoallv, 1, MPI_INT, MPI_COMM_WORLD);

    // calculate the place
//    int* send_displs = (int*)calloc(buckets_number, sizeof(int));
//    int* recv_displs = (int*)calloc(buckets_number, sizeof(int));
//    // send_diapls[0] and recv_displs[0] are both equal to 0, so here
//    send_displs[0] = 0;
//    recv_displs[0] = 0;
//    for (int i = 1; i < buckets_number; i++){
//        send_displs[i] = i * curr_proc_data_size;
//        recv_displs[i] = recv_displs[i-1] + recv_count_alltoallv[i-1];
//    }
//
//    cout << "Line: " << __LINE__ << ", Before alltoallv send_displs and recv_displs, the rank of this processor is " << curr_rank << endl;
//    display_int_array(send_displs, buckets_number);
//    display_int_array(recv_displs, buckets_number);
//
//    // use alltoallv to communicate numbers in each processores
//    recv_bucket_alltoallv = (unsigned long*)calloc(number_size, sizeof(unsigned long));
//    MPI_Alltoallv(bucket, nitems, send_displs, MPI_LONG, recv_bucket_alltoallv, recv_count_alltoallv, recv_displs, MPI_LONG, MPI_COMM_WORLD);
//
//    cout << "Line: " << __LINE__ << ", After alltoallv, the rank of this processor is " << curr_rank << endl;
//    display(recv_bucket_alltoallv, number_size);
//    display_int_array(recv_displs, buckets_number);
//
//    // step 5, each process sorts its own numbers.
//    unsigned long *result; // receive all numbers per process
//    int recv_total_count = 0; // total number of each process
//    for (int i = 0; i < buckets_number; i ++)
//        recv_total_count += recv_count_alltoallv[i];
//
//    result = new unsigned long[recv_total_count];
//    int count = 0;
//    for(unsigned long i=0;i<number_size;i++)
//    {
//        if(recv_bucket_alltoallv[i]!=INF)
//        {
//            result[count++] = recv_bucket_alltoallv[i];
//        }
//    }
//
//    cout << "Line: " << __LINE__ << ", display result, rank : " << curr_rank << endl;
//    display(result, recv_total_count);
//
//    // just use qsort of stdlib
//    qsort(result, recv_total_count, sizeof(unsigned long), IncOrder);
//
//    cout << "Line: " << __LINE__ << ", display sorted result, rank : " << curr_rank << endl;
//    display(result, recv_total_count);
//
//    //step 6, Gather the results to rank 0
//    int *recv_cnt = new int[buckets_number]; // receive count from processes
//    unsigned long *sorted = new unsigned long[number_size]; // sorted numbers
//    int *final_displs = new int[buckets_number]; // final displaces
//
//    // gather the count from each process.
//    MPI_Gather(&recv_total_count, 1, MPI_INT, recv_cnt, 1, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);
//    final_displs[0]=0;
//    for(int i = 1; i < buckets_number; i++)
//    {
//        final_displs[i] = recv_displs[i-1] + recv_cnt[i-1];
//    }
//
//    cout << "Line: " << __LINE__ << ", display final_displs, rank : " << curr_rank << endl;
//    display_int_array(final_displs, processes_number);
//
//    MPI_Gatherv(result, recv_total_count, MPI_LONG, sorted, recv_cnt, final_displs, MPI_LONG, MASTER_RANK, MPI_COMM_WORLD);
//    cost_time += MPI_Wtime();
//    cout << "Line: " << __LINE__ << ", time of curr_rank " << curr_rank << " : " << cost_time<<endl;
//
//    //print the sorted data on curr_rank 0
//	if(curr_rank == MASTER_RANK)
//	{
//        cout << "Line: " << __LINE__ << ", display sorted result, rank : " << curr_rank << endl;
//        display(sorted, number_size);
//	}

    MPI_Finalize();
    return 0;
}