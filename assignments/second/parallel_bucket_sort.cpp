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
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <time.h>
#include <cstdlib>
#include <math.h>

// ** note ** in my machine cannot find tr1, but in mighty there is a runtime error if no tr1

#define INF (-1)
using namespace std;

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


float *generate_random_number (int xmin, int xmax, int size){

    float *data = (float *)malloc(size * sizeof(float));
    for (int i = 0; i < size; i++)
        data[i] = drand48() * (xmax - xmin - 1) + xmin;

    return data;
}

// comparator for float
int float_comparator(const void *x1, const void *x2) {
    float *f1 = (float *)x1;
    float *f2 = (float *)x2;
    float diff = *f1 - *f2;

    return (diff < 0) ? -1 : 1;

//    return (*((float *)x1)-*((float *)x2)) < 0 ? -1 : 1; //(diff < 0) ? -1 : 1;
}

//IncOrder for qsort
int IncOrder(const void *e1, const void *e2)
{
    return (*((unsigned long *)e1)-*((unsigned long *)e2));
}

void display(float *array, float size) {
    for(int i = 0; i<size; i++)
        cout << array[i] << " ";
    cout << endl;
}

void display_int_array_with_info(int *array, int size, int curr_rank, string info) {
    cout << "curr_rank = " << curr_rank << ", " << info ;
    for(int i = 0; i<size; i++)
        cout << array[i] << " ";
    cout << endl;
}

void display_int_array(int *array, int size, int curr_rank) {
    cout << "curr_rank = " << curr_rank << endl;
    for(int i = 0; i<size; i++)
        cout << array[i] << " ";
    cout << endl;
}

/*
void display_int_array(int *array, int size, int curr_rank, char info) {
    cout << "curr_rank = " << curr_rank << endl;
    for(int i = 0; i<size; i++)
        cout << array[i] << " ";
    cout << endl;
}*/

int main(int argc, char **argv)
{
    if(argc!=2)
    {
        cout<<"Please check your input"<<endl;
        exit(0);
    }

    double cost_time = - MPI_Wtime();

    const float xmin = 10.0;
    const float xmax = 100; //250000;
    const int MASTER_RANK = 0; // the master's rank
    const bool IS_DEBUG = true;

    //rank of the current process, the number of processes, size of data on the current process, the number of numbers
    int curr_rank, processes_number, curr_proc_data_size, number_size, i;

    //data on the current process, the unsorted data, the buckets in the current process, the bucket after alltoallv function
    float *curr_proc_data, *original_data, *small_buckets, *recv_big_bucket;

    //Blocks until all processes in the communicator have reached this routine
    //    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Init(&argc, &argv); //initial
    MPI_Comm_rank(MPI_COMM_WORLD, &curr_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processes_number);
    number_size = atoll(argv[1]); // the size of numbers need to be sorted, an input parameter
    curr_proc_data_size = number_size / processes_number; // the number of numbers is total_number / processers

    if (IS_DEBUG)
        cout << "number_size = " << number_size << ", curr_proc_data_size = " << curr_proc_data_size << endl;

    // step 1, initial the number
    if(curr_rank == MASTER_RANK)
    {
        original_data = generate_random_number (xmin, xmax, number_size);
    }

    // step 2, scatter evenly to all processes.
    curr_proc_data = (float *)malloc(curr_proc_data_size * sizeof(float));//new float[curr_proc_data_size];
//    final_sorted_data=new long[number_size];
    MPI_Scatter(original_data, curr_proc_data_size, MPI_FLOAT, curr_proc_data, curr_proc_data_size, MPI_FLOAT, MASTER_RANK, MPI_COMM_WORLD);

    if (IS_DEBUG) {
        for (i = 0; i < curr_proc_data_size; i++) {
            cout << "after scatter rank " << curr_rank << " : curr_proc_data [i] = " << curr_proc_data[i] << cost_time
                 << endl;
        }
    }
    // OK

    // step 3, each processor loops each number to determine which small bucket they should in.
    //initialize small_buckets
    small_buckets = (float *)malloc(number_size * sizeof(float));//new float[number_size];
    for(i = 0; i < number_size; i++)
    {
        small_buckets[i] = INF;
    }

    // update small_buckets
    // curr_proc_data --> small_buckets
    int buckets_number = processes_number;

    //initialize number of items, used to storte the size of numbers in small buckets
    int *nitems = (int *)malloc(buckets_number * sizeof(int));
    for (int i = 0; i < buckets_number; ++i)
        nitems[i] = 0;
    float step = (xmax - xmin)/buckets_number;
    if (IS_DEBUG)
        cout << __LINE__ << ", curr_rank = " << curr_rank << ", "  << ", step = " << step << endl;

    small_buckets = (float*)calloc(buckets_number * curr_proc_data_size, sizeof(float));
    for(i = 0; i < number_size; i++)
        small_buckets[i] = INF;

//    float *big_bucket = calloc(ntotal, sizeof(float *));
//    for (i = 0; i < ntotal; ++i) bucket[i] = 0;

    for (i = 0; i < curr_proc_data_size; i++)
    {
        int bktno = (int)floor((curr_proc_data[i] - xmin)/step);// in which small bucket
        int idx = bktno * curr_proc_data_size + nitems[bktno];// index in the bucket
        if (IS_DEBUG)
            cout << __LINE__ << "curr_rank = " << curr_rank << ", curr_proc_data[i] = " << curr_proc_data[i] << ", bktno = " << bktno << ", idx = " << idx << endl;
        small_buckets[idx] = curr_proc_data[i];
        ++nitems[bktno];
    }

    if (IS_DEBUG) {
        cout << __LINE__ << ", curr_rank = " << curr_rank << ", " << endl;
        display_int_array_with_info(nitems, buckets_number, curr_rank, "nitems");
    }

    // step 4, each processor scatter its numbers to proper processors and gather its own proper numbers from others
    // firstly, need to let all processores know how many numbers should recv from each processor
    int* recv_count_alltoallv = (int*)calloc(buckets_number, sizeof(int));
 //   int send_count, recv_count = 1;
    MPI_Alltoall(nitems, 1, MPI_INT, recv_count_alltoallv, 1, MPI_INT, MPI_COMM_WORLD);

    if (IS_DEBUG){
    //    cout << __LINE__ << ", curr_rank = " << curr_rank << ", " << endl;
        display_int_array(recv_count_alltoallv, buckets_number, curr_rank);
    }
    // calculate the place
    int* send_displs = (int*)calloc(buckets_number, sizeof(int));
    int* recv_displs = (int*)calloc(buckets_number, sizeof(int));
    // send_diapls[0] and recv_displs[0] are both equal to 0, so here
    send_displs[0] = 0;
    recv_displs[0] = 0;
    for (i = 1; i < buckets_number; i++){
        send_displs[i] = i * curr_proc_data_size;
        recv_displs[i] = recv_displs[i-1] + recv_count_alltoallv[i-1];
    }

    if (IS_DEBUG) {
        cout << "Line: " << __LINE__ << ", Before alltoallv send_displs and recv_displs, the rank of this processor is "
             << curr_rank << endl;
        display_int_array(send_displs, buckets_number, curr_rank);
        display_int_array(recv_displs, buckets_number, curr_rank);
    }

    // use alltoallv to communicate numbers in each processores
    recv_big_bucket = (float *)calloc(number_size, sizeof(float));
    MPI_Alltoallv(small_buckets, nitems, send_displs, MPI_FLOAT, recv_big_bucket, recv_count_alltoallv, recv_displs, MPI_FLOAT, MPI_COMM_WORLD);

    if (IS_DEBUG) {
        cout << "Line: " << __LINE__ << ", After alltoallv, the rank of this processor is " << curr_rank << endl;
        display(recv_big_bucket, number_size);
        display_int_array(recv_displs, buckets_number, curr_rank);
    }

    // step 5, each process sorts its own numbers.
    float *result; // receive all numbers per process
    int recv_total_count = 0; // total number of each process
    for (i = 0; i < buckets_number; i ++)
        recv_total_count += recv_count_alltoallv[i];

    result = (float *)calloc(recv_total_count, sizeof(float));//new float[recv_total_count];
    int count = 0;
    for(i = 0; i < number_size; i++)
    {
        if(recv_big_bucket[i]!=INF)
            result[count++] = recv_big_bucket[i];
    }

    if (IS_DEBUG) {
        cout << "Line: " << __LINE__ << ", display result, rank : " << curr_rank << endl;
        display(result, recv_total_count);
    }

    // just use qsort of stdlib
    qsort(result, recv_total_count, sizeof(float), float_comparator);

    if (IS_DEBUG) {
        cout << "Line: " << __LINE__ << ", display sorted result, rank : " << curr_rank << endl;
        display(result, recv_total_count);
    }

    //step 6, Gather the results to rank 0
//    int *recv_cnt = new int[buckets_number]; // receive count from processes
    float *sorted = (float*)calloc(number_size, sizeof(float)); //new float[number_size]; // sorted numbers
    int *final_displs = (int*)calloc(buckets_number, sizeof(int));//new int[buckets_number]; // final displaces

    if (IS_DEBUG) {
        cout << "Line: " << __LINE__ << ", recv_total_count : " << recv_total_count << endl;
    }

    // gather the count from each process.
//    memset(recv_count_alltoallv, 0, buckets_number*sizeof(int));
    MPI_Gather(&recv_total_count, 1, MPI_INT, recv_count_alltoallv, 1, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);

    if (IS_DEBUG) {
        cout << "Line: " << __LINE__ << ", gather count result, rank : " << curr_rank << endl;
        display_int_array(recv_count_alltoallv, buckets_number, curr_rank);
    }
    // it is ok

    final_displs[0]=0;
    for(i = 1; i < buckets_number; i++)
    {// recv_displs --> recv_count_alltoallv
        final_displs[i] = final_displs[i-1] + recv_count_alltoallv[i-1]; //recv_cnt[i-1];
    }

    if (IS_DEBUG) {
        cout << "Line: " << __LINE__ << ", display final_displs, rank : " << curr_rank << endl;
        display_int_array(final_displs, processes_number, curr_rank);
    }

    //result ? big_bucket ?
    // recv_cnt --> recv_count_alltoallv
    MPI_Gatherv(result, recv_total_count, MPI_FLOAT, sorted, recv_count_alltoallv, final_displs, MPI_FLOAT, MASTER_RANK, MPI_COMM_WORLD);
    cost_time += MPI_Wtime();
    cout << "Line: " << __LINE__ << ", time of curr_rank " << curr_rank << " : " << cost_time<<endl;

    //print the sorted data on curr_rank 0
	if(IS_DEBUG && curr_rank == MASTER_RANK)
	{
        cout << "Line: " << __LINE__ << ", display sorted result, rank : " << curr_rank << endl;
        display(sorted, number_size);
	}

    MPI_Finalize();
    return 0;
}