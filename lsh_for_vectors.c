#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <limits.h>
#include "lsh_for_vectors.h"
#include "lsh_funcs.h"


double** grid_to_vector(double** curves, double delta, int input_items_counter, int dimension, float* t)
{
    double** vector = malloc(sizeof(double*) * input_items_counter);
    for(int i=0; i<input_items_counter; i++)
    {
        vector[i] = malloc(sizeof(double) * dimension);
        for(int j=0; j<dimension; j++)
            vector[i][j] = 3*dimension;
    }

    double G[dimension];
    int flag, count;
    for(int i=0; i<input_items_counter; i++)
    {
        count = 0;
        for(int j=0; j<dimension; j++)
        {
            flag = 0;
            G[j] = floor((curves[i][j] / delta) + 1/2) * delta;     // operation for grid
            G[j] = G[j] + t[j];    // shift vector by t
            
            if(count > 0)
            {
                if((int)vector[i][count-1] == (int)G[j])
                {
                    flag = 1;
                }
            }
            if(flag == 0)
            {
                vector[i][count] = G[j];
                count++;
            }
        }
    }

    return vector;  // vector has the new vectors that are going to be in lsh
}

// void lsh_for_vectors(double** vectors, int input_items_counter, int dimension, int k, int L, char** names, FILE* query_file_ptr, FILE* output_file_ptr, double delta, double** curves, float* t, int n)
void lsh_for_vectors(double** vectors, int input_items_counter, char** names, double** query_vectors, int query_items_counter, char** query_names, int dimension, double** curves, double** query_curves, FILE* output_file_ptr, int k, double delta, float* t, int n)
{
    float** h_p_result = malloc(sizeof(float*) * input_items_counter); // array with the results of the h function
    for(int i=0; i<input_items_counter; i++)
        h_p_result[i] = malloc(sizeof(float) * k);
    for(int i=0; i<input_items_counter; i++)
    {
        srand(time(0));
        for(int j=0; j<k; j++)
        {
            h_p_result[i][j] = h_function(vectors, i, dimension);   // h_function
        }
    }

    // Hash table for input file
    int hash_index;
    int TableSize = input_items_counter / 8;
    int M_hash = (int)pow(2, 32) - 5;
    struct Hash_Node* hash_tables[TableSize];
    for(int i=0; i<TableSize; i++)
    {
        hash_tables[i] = NULL;
    }

    int r[k];
    int ID;
    for(int i=0; i<k; i++)
    {
        r[i] = rand() % 10;  // r for g function
    }
    for(int i=0; i<input_items_counter; i++)
    {
        hash_index = 0;
        for(int j=0; j<k; j++)
        {
            hash_index = hash_index + (int)h_p_result[i][j] * r[j];  // g function
        }
        hash_index = hash_index % M_hash;    // mod M
        if(hash_index < 0)
        {
            hash_index = hash_index * (-1);     // only positive number
        }
        ID = hash_index;
        hash_index = hash_index % TableSize;    // mod TableSize
        struct Hash_Node* data_item = (struct Hash_Node*)malloc(sizeof(struct Hash_Node));  // new node
        data_item->name = malloc(sizeof(char*) + 1);
        data_item->name = names[i];    // fill it
        data_item->item = i;
        data_item->ID = ID;
        if(hash_tables[hash_index] == NULL)  // put it in the list if list is empty
        {
            hash_tables[hash_index] = data_item;
            hash_tables[hash_index]->ID = ID;
        }
        else    // put it after the last node of the list
        {
            struct Hash_Node* temp = hash_tables[hash_index];
            while(hash_tables[hash_index]->next != NULL)
            {
                hash_tables[hash_index] = hash_tables[hash_index]->next;
            }
            hash_tables[hash_index]->next = data_item;
            hash_tables[hash_index]->next->ID = ID;
            hash_tables[hash_index] = temp;
        }
    }

    float** h_q_result = malloc(sizeof(float*) * query_items_counter); // array with the results of the h function
    for(int i=0; i<query_items_counter; i++)
        h_q_result[i] = malloc(sizeof(float) * k);
    for(int i=0; i<query_items_counter; i++)
    {
        srand(time(0));
        for(int j=0; j<k; j++)
        {
            h_q_result[i][j] = h_function(query_vectors, i, dimension);   // h_function
        }
    }

    // find the nearest neighbor of each query
    char* nearest_neighbor = malloc(sizeof(char*) + 1);
    fprintf(output_file_ptr, "Hash Table no%d\n", n);
    for(int m=0; m<query_items_counter; m++)    // for every query show the results
    {
        hash_index = 0;
        for(int j=0; j<k; j++)
        {
            hash_index = hash_index + (int)h_q_result[m][j] * r[j]; // find the bucket of the query (g function)
        }
        hash_index = hash_index % M_hash;    // mod M
        if(hash_index < 0)
        {
            hash_index = hash_index * (-1);
        }
        int k_ID = hash_index;  // ID for comparison
        hash_index = hash_index % TableSize;    // mod TableSize

        float min_dist = 1000000.0;
        // struct timeval start, stop;
        // gettimeofday(&start, 0);
        while(hash_tables[hash_index] != NULL)
        {
            if(k_ID == hash_tables[hash_index]->ID)  // compare the IDs
            {
                // calculate distance
                float dist = 0;
                for(int d=0; d<dimension; d++)
                {
                    dist = dist + pow((query_curves[m][d] - curves[hash_tables[hash_index]->item][d]), 2);   // distance of items
                }
                dist = sqrt(dist);
                if(dist < min_dist) // minimun LSH distance
                {
                    min_dist = dist;
                    nearest_neighbor = names[hash_tables[hash_index]->item];
                }
            }
            hash_tables[hash_index] = hash_tables[hash_index]->next;        
        }
        // gettimeofday(&stop, 0);
        // long sec = stop.tv_sec - start.tv_sec;
        // long mic_sec = stop.tv_usec - start.tv_usec;
        // double lsh_time = sec + mic_sec*1e-6;

        // calculate true distance and time
        float true_min_dist = 1000000.0;
        char* true_nearest_neighbor = malloc(sizeof(char*) + 1);
        // gettimeofday(&start, 0);
        for(int i=0; i<input_items_counter; i++)
        {
            float dist = 0;
            for(int d=0; d<dimension; d++)
            {
                dist = dist + pow((query_curves[m][d] - curves[i][d]), 2);
            }
            dist = sqrt(dist);
            if(dist < true_min_dist && dist >= 0) // minimun LSH distance
            {
                true_min_dist = dist;
                true_nearest_neighbor = names[i];
            }
        }
        // gettimeofday(&stop, 0);
        // sec = stop.tv_sec - start.tv_sec;
        // mic_sec = stop.tv_usec - start.tv_usec;
        // double true_time = sec + mic_sec*1e-6;

        if(min_dist != 1000000.0)
        {
            // print query
            fprintf(output_file_ptr, "Query: %s\n", query_names[m]);

            // print algorithm
            fprintf(output_file_ptr, "Algorithm: LSH_Vector\n");

            // print nearest neighbor
            fprintf(output_file_ptr, "Approximate Nearest neighbor: %s\n", nearest_neighbor);

            // print true nearest neighbor
            fprintf(output_file_ptr, "True Nearest neighbor: %s\n", true_nearest_neighbor);

            // print LSH distance
            fprintf(output_file_ptr, "distanceApproximate: %f\n", min_dist);

            // print true distance
            fprintf(output_file_ptr, "distanceTrue: %f\n", true_min_dist);
        }
        else
        {
            fprintf(output_file_ptr, "Query: %s\n", query_names[m]);
            fprintf(output_file_ptr, "There is not a near vector in this bucket!\n");
        }
        fprintf(output_file_ptr, "\n");
    }
    fprintf(output_file_ptr, "AAAAAAAAAA\n");
    // }
}