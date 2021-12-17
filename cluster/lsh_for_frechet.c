#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <limits.h>
#include "lsh_for_frechet.h"
#include "lsh_funcs.h"

float distance_computation(float** distance, int dimension, double** curves, double** query_curves, int m, int item)
{
    distance[0][0] = sqrt(pow(curves[item][0] - query_curves[m][0], 2));

    for(int i=1; i<dimension; i++)
    {
        distance[i][0] = fmax(distance[i-1][0], sqrt(pow(curves[item][i] - query_curves[m][0], 2)));
    }
    for(int i=1; i<dimension; i++)
    {
        distance[0][i] = fmax(distance[0][i-1], sqrt(pow(curves[item][0] - query_curves[m][i], 2)));
    }

    for(int i=1; i<dimension; i++)
    {
        for(int j=1; j<dimension; j++)
        {
            distance[i][j] = fmax(fmin(fmin(distance[i-1][j], distance[i][j-1]), distance[i-1][j-1]), sqrt(pow(curves[item][i] - query_curves[m][j], 2)));
        }
    }

    return distance[dimension-1][dimension-1];
}

double** grid_to_frechet(double** curves, double delta, int input_items_counter, int dimension, float* t)
{
    struct Node* grid_array[input_items_counter][dimension];
    for(int i=0; i<input_items_counter; i++)
    {
        for(int j=0; j<dimension; j++)
        {
            grid_array[i][j] = malloc(sizeof(struct Node));
            grid_array[i][j]->x = 3*dimension;  // for padding
            grid_array[i][j]->y = 3*dimension;  // fill it and change only the different ones
        }
    }
    int flag, count;
    double G[dimension];
    int G_index[dimension];
    for(int i=0; i<input_items_counter; i++)
    {
        count = 0;
        for(int j=0; j<dimension; j++)
        {
            flag = 0;
            G[j] = floor((curves[i][j] + 1/2) / delta) * delta;     // snapping
            G[j] = G[j] + t[j];    // shift vector by t
            G_index[j] = j * delta;

            // if something is the same then throw it
            if(count > 0)
            {
                if((int)grid_array[i][count-1]->y == (int)G[j])
                    flag = 1;
            }

            if(flag == 0)   // else fill the array with the results
            {
                grid_array[i][count]->x = G_index[j];
                grid_array[i][count]->y = G[j];
                count++;
            }
        }

        // if(i==0)
        // {
        //     for(int p=0; p<dimension; p++)
        //         printf("%d %f\n", grid_array[i][p]->x, grid_array[i][p]->y);
        // }
    }

    int new_dimension = 2 * dimension;
    double** vector = malloc(sizeof(double*) * input_items_counter);
    for(int i=0; i<input_items_counter; i++)
        vector[i] = malloc(sizeof(double) * new_dimension);

    for(int i=0; i<input_items_counter; i++)
    {
        count = 0;
        for(int j=0; j<dimension; j++)
        {
            vector[i][count] = grid_array[i][j]->x;  // vector is the array that has every clue of the struct aligned
            count++;
            vector[i][count] = grid_array[i][j]->y;
            count++;
        }
    }

    // for(int i=0; i<new_dimension; i++)
    //     printf("%f\n", vector[0][i]);

    return vector;  // vector has the new vectors that are going to be in lsh
}

void lsh_for_frechet(double** vectors, int input_items_counter, char** names, double** query_vectors, int query_items_counter, char** query_names, int dimension, double** curves, double** query_curves, FILE* output_file_ptr, int k, int n)
{
    float** h_p_result = malloc(sizeof(float*) * input_items_counter); // array with the results of the h function
    for(int i=0; i<input_items_counter; i++)
        h_p_result[i] = malloc(sizeof(float) * k);
    for(int i=0; i<input_items_counter; i++)
    {
        srand(time(0));
        for(int j=0; j<k; j++)
        {
            h_p_result[i][j] = h_function(vectors, i, 2*dimension);   // h_function
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
        }
        else    // put it after the last node of the list
        {
            struct Hash_Node* temp = hash_tables[hash_index];
            while(hash_tables[hash_index]->next != NULL)
            {
                hash_tables[hash_index] = hash_tables[hash_index]->next;
            }
            hash_tables[hash_index]->next = data_item;
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
            h_q_result[i][j] = h_function(query_vectors, i, 2*dimension);   // h_function
        }
    }

    float** distance = malloc(sizeof(float*) * dimension);
    for(int i=0; i<dimension; i++)
        distance[i] = malloc(sizeof(float) * dimension);


    char* nearest_neighbor = malloc(sizeof(char*) + 1);
    float dist;
    // find the nearest neighbor of each query
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
        dist = min_dist;
        while(hash_tables[hash_index] != NULL)
        {
            if(k_ID == hash_tables[hash_index]->ID)  // compare the IDs
            {
                // calculate distance
                for(int i=0; i<dimension; i++)
                    for(int j=0; j<dimension; j++)
                        distance[i][j] = -1;
                dist = distance_computation(distance, dimension, curves, query_curves, m, hash_tables[hash_index]->item);
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
        dist = true_min_dist;
        // gettimeofday(&start, 0);
        for(int i=0; i<input_items_counter; i++)
        {
            for(int i=0; i<dimension; i++)
                for(int j=0; j<dimension; j++)
                    distance[i][j] = -1;
            dist = distance_computation(distance, dimension, curves, query_curves, m, i);
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
}