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

unsigned int number_size;//length of the unsorted data
int curr_rank;//rank of the current process
int processes_number;//number of processes
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
    return (*((long *)e1)-*((long *)e2));
}

int main(int argc, char **argv)
{
    if(argc!=2)
    {
        cout<<"Please check your input"<<endl;
        exit(0);
    }
    long curr_proc_data_size;//size of data on the current process
    long *proc_data;//data on the current process
    long *original;//the unsorted data
//    long *final_sorted_data;//the sorted data
    long *pivot_list;//the pivot list
    long *proc_buckets;//the buckets in the current process
    long *final_buckets;//the bucket after alltoall function
    double cost_time;

    MPI_Init(&argc, &argv); //initial

    // TODO  is it need?
    MPI_Barrier(MPI_COMM_WORLD); //Blocks until all processes in the communicator have reached this routine
    cost_time=-MPI_Wtime();
    MPI_Comm_rank(MPI_COMM_WORLD, &curr_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processes_number);

    if(0==curr_rank)
    {
        number_size = atoll(argv[1]);

        original = generate_array_with_random (number_size, number_size);
        curr_proc_data_size = number_size/processes_number;
        cout << curr_proc_data_size << endl;
    }
    MPI_Bcast(&number_size, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&curr_proc_data_size, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    proc_data=new long[curr_proc_data_size];

//    final_sorted_data=new long[number_size];
    pivot_list=new long[processes_number];
    MPI_Scatter(original, curr_proc_data_size, MPI_LONG, proc_data, curr_proc_data_size, MPI_LONG, 0, MPI_COMM_WORLD);

    //initial local sort
    qsort(proc_data, curr_proc_data_size, sizeof(long), IncOrder);

    MPI_Gather(proc_data, 1, MPI_LONG, pivot_list, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    if(0==curr_rank)
    {
        qsort(pivot_list, processes_number, sizeof(long), IncOrder);
    }

    proc_buckets=new long[number_size];

    //initialize proc_buckets
    for(long i=0;i<number_size;i++)
    {
        proc_buckets[i]=INF;
    }

    long *index=new long[processes_number];

    //initialize index
    for(int i=0;i<processes_number;i++)
    {
        index[i]=0;
    }
    MPI_Bcast(pivot_list, processes_number, MPI_LONG, 0, MPI_COMM_WORLD);


    //update proc_buckets
    for(long i=0;i<curr_proc_data_size;i++)
    {
        for(int j=0;j<processes_number-1;j++)
        {
            if(proc_data[i]>=pivot_list[j]&&proc_data[i]<pivot_list[j+1])
            {
                proc_buckets[j*curr_proc_data_size+index[j]]=proc_data[i];
                index[j]=index[j]+1;
            }
        }
        if(proc_data[i]>=pivot_list[processes_number-1])
        {
            proc_buckets[(processes_number-1)*curr_proc_data_size+index[processes_number-1]]=proc_data[i];
            index[processes_number-1]=index[processes_number-1]+1;
        }
    }

    //creation of a new datatype BUCKETS
    MPI_Datatype BUCKETS;
    MPI_Type_contiguous(curr_proc_data_size, MPI_LONG, &BUCKETS);
    MPI_Type_commit(&BUCKETS);


    final_buckets=new long[processes_number*curr_proc_data_size];
    for(long i=0;i<number_size;i++)
    {
        final_buckets[i]=0;
    }

    //the alltoall function to get the final buckets in processes
    MPI_Alltoall(proc_buckets, 1, BUCKETS, final_buckets, 1, BUCKETS, MPI_COMM_WORLD);

    MPI_Type_free(&BUCKETS);

    long *result;
    long count=0;
    for(long i=0;i<number_size;i++)
    {
        if(final_buckets[i]!=INF)
        {
            count++;
        }
    }
    result=new long[count];
    count=0;
    for(long i=0;i<number_size;i++)
    {
        if(final_buckets[i]!=INF)
        {
            result[count++]=final_buckets[i];
        }
    }

    qsort(result, count, sizeof(long), IncOrder);


    //Gather the results to rank 0
    int *recv_cnt=new int[processes_number];
    long *sorted=new long[number_size];
    int *displs=new int[processes_number];

    MPI_Gather(&count, 1, MPI_LONG, recv_cnt, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    displs[0]=0;
    for(int i=1;i<processes_number;i++)
    {
        displs[i]=displs[i-1]+recv_cnt[i-1];
    }

    MPI_Gatherv(result, count, MPI_LONG, sorted, recv_cnt, displs, MPI_LONG, 0, MPI_COMM_WORLD);
    cost_time+=MPI_Wtime();
    cout<<"time of curr_rank "<<curr_rank<<":"<<cost_time<<endl;

    //print the sorted data on curr_rank 0
	if(0==curr_rank)
	{
		for(long i=0;i<number_size;i++)
		{
			cout<<sorted[i]<<" ";
		}
		cout<<endl;
	}

    MPI_Finalize();
    return 0;
}