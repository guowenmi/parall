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
    unsigned long *final_buckets;//the bucket after alltoallv function
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

    //initialize index, used to record the size of numbers in small buckets
    unsigned long *index = new unsigned long[processes_number];
    for(int i = 0; i < processes_number; i++)
    {
        index[i] = 0;
    }

    // update small_buckets
    // curr_proc_data --> small_buckets
    int buckets_number = processes_number;
    unsigned long *bucket = calloc(buckets_number * curr_proc_data_size, sizeof(long));

    //initialize number of items, used to storte the size of numbers in small buckets
    int *nitems = calloc(buckets_number, sizeof(int));
    unsigned long step = number_size/processes_number;

    for (int i = 0; i < curr_proc_data_size; i++)
    {
        int bktno = floor(curr_proc_data[i]/step);// in which bucket
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
    int* send_displs = (int*)calloc(buckets_number, sizeof(int));
    int* recv_displs = (int*)calloc(buckets_number, sizeof(int));
    for (int i = 1; i < buckets_number; i++){
        send_displs[i] = i * curr_proc_data;
        recv_displs[i] = recv_displs[i-1]+recv_count_alltoallv[i-1];
    }

    // use alltoallv to communicate numbers in each processores
    unsigned long* big_bucket = calloc(number_size, sizeof(unsigned long));
    MPI_Alltoallv(bucket, nitems, send_displs, MPI_LONG, big_bucket, recv_count_alltoallv, recv_displs, MPI_LONG, MPI_COMM_WORLD);

    cout << "the rank of this processor is " << curr_rank << endl;
    display(big_bucket, recv_count_alltoallv);

//
//    vector<unsigned long> bucket[processes_number];
//    int position = 0;
//    for(unsigned long i = 0; i<curr_proc_data_size; i++)  {          //put elements into different buckets
//        position = floor (curr_proc_data[i] / number_size * (processes_number - 1));
//        bucket[position].push_back(curr_proc_data[i]);
//    }
//
//    unsigned long index = 0;
//    for(int i = 0; i < processes_number; i++) {
//        while(!bucket[i].empty()) {
//            curr_proc_data[index++] = *(bucket[i].begin());
//            bucket[i].erase(bucket[i].begin());
//        }
//    }
//
//    display(arr, n);


//
//    // from bucket.c
//    int i, count;
//
//    // The range covered by one bucket
//    float stepsize = number_size / processes_number ; //(x2 - x1) / nbuckets;
//
//    // The number of items thrown into each bucket. We would expect each
//    // bucket to have a similar number of items, but they won't be
//    // exactly the same. So we keep track of their numbers here.
//    int *nitems = malloc(processes_number * sizeof(int));
//    for (i = 0; i < processes_number; ++i) nitems[i] = 0;
//
//    // Toss the data items into the correct bucket
//    for (i = 0; i < curr_proc_data_size; ++i) {
//
//        // What bucket does this data value belong to?
//        int bktno = (int) floor(curr_proc_data[i]  / stepsize);
//        int idx = bktno * curr_proc_data_size + nitems[bktno];
//
//        // Put the data value into this bucket
//        bucket[idx] = data[i];
//        ++nitems[bktno];
//    }
//    // from bucket.c


//    MPI_Bcast(pivot_list, processes_number, MPI_LONG, 0, MPI_COMM_WORLD);

    //update small_buckets
//    for(unsigned long i = 0; i < curr_proc_data_size; i++)
//    {
//        for(int j = 0; j < processes_number - 1; j++)
//        {
//            if(curr_proc_data[i] >= pivot_list[j] && curr_proc_data[i] < pivot_list[j+1])
//            {
//                small_buckets[j*curr_proc_data_size+index[j]] = curr_proc_data[i];
//                index[j] = index[j] + 1;
//            }
//        }
//
//        if(curr_proc_data[i]>=pivot_list[processes_number-1])
//        {
//            small_buckets[(processes_number-1)*curr_proc_data_size+index[processes_number-1]]=curr_proc_data[i];
//            index[processes_number-1]=index[processes_number-1]+1;
//        }
//    }
//
//    //creation of a new datatype BUCKETS
//    MPI_Datatype BUCKETS;
//    MPI_Type_contiguous(curr_proc_data_size, MPI_LONG, &BUCKETS);
//    MPI_Type_commit(&BUCKETS);
//
//    final_buckets=new unsigned long[processes_number*curr_proc_data_size];
//    for(unsigned long i=0;i<number_size;i++)
//    {
//        final_buckets[i]=0;
//    }
//
//    //the alltoall function to get the final buckets in processes
//    MPI_Alltoall(small_buckets, 1, BUCKETS, final_buckets, 1, BUCKETS, MPI_COMM_WORLD);
//    MPI_Type_free(&BUCKETS);


    // step 5, each process sorts its own numbers.
    unsigned long *result;
    unsigned long count=0;
    for(unsigned long i=0;i<number_size;i++)
    {
        if(final_buckets[i]!=INF)
        {
            count++;
        }
    }

    result = new unsigned long[count];
    count = 0;
    for(unsigned long i=0;i<number_size;i++)
    {
        if(final_buckets[i]!=INF)
        {
            result[count++] = final_buckets[i];
        }
    }

    qsort(result, count, sizeof(unsigned long), IncOrder);


    //step 6, Gather the results to rank 0
    int *recv_cnt = new int[processes_number];
    unsigned long *sorted = new unsigned long[number_size];
    int *displs = new int[processes_number];

    MPI_Gather(&count, 1, MPI_LONG, recv_cnt, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    displs[0]=0;
    for(int i=1;i<processes_number;i++)
    {
        displs[i]=displs[i-1]+recv_cnt[i-1];
    }

    MPI_Gatherv(result, count, MPI_LONG, sorted, recv_cnt, displs, MPI_LONG, 0, MPI_COMM_WORLD);
    cost_time += MPI_Wtime();
    cout << "time of curr_rank " << curr_rank << " : " << cost_time<<endl;

    //print the sorted data on curr_rank 0
	if(curr_rank == MASTER_RANK)
	{
		for(long i = 0; i < number_size; i++)
		{
			cout << sorted[i] << " ";
		}
		cout << endl;
	}

    MPI_Finalize();
    return 0;
}