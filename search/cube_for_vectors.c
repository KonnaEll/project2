#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <limits.h>
#include "cube_for_vectors.h"
#include "cube_funcs.h"
#include "lsh_funcs.h"


void cube_for_vectors(int input_items_counter, char** names, int query_items_counter, char** query_names, int dimension, double** curves, double** query_curves, FILE* output_file_ptr, int k, int probes,int M)
{
    // create 2 tables with 2 dimension  and we put inside 0,1//
    int** f = (int**)malloc(input_items_counter * sizeof(int*));
    for(int i=0; i<input_items_counter; i++)
    {
        f[i] = (int*)malloc(sizeof(int) * k);
    }
    for (int i=0;i<input_items_counter; i++){
        for (int j=0;j<k;j++){
            f[i][j]=rand() % 2;
        }
    }


    float** h_p_result = malloc(sizeof(float*) * input_items_counter * sizeof(int*)); // array with the results of the h function
    for(int i=0; i<input_items_counter; i++)
        h_p_result[i] = malloc(sizeof(float) * k);
    for(int i=0; i<input_items_counter; i++)
    {
        srand(time(0));
        for(int j=0; j<k; j++)
        {
            h_p_result[i][j] = h_function(curves, i, dimension);   // h_function
            for(int a=0;a<i;a++){
                for(int b=0;b<j;b++){
                    if(h_p_result[i][j] == h_p_result[a][b]){   //if we find the same result we put the number 0 in f table//
                        f[i][j]=0;
                    }
                }
            }
        }
    }
    int * a=h_p(f,k,input_items_counter);

    // Hash table for input file
    int TableSize = pow(2,k);
    int hash_index;
    struct Hash_Node* hash_tables[TableSize];
    for(int i=0; i<TableSize; i++)
    {
        hash_tables[i] = NULL;
    }

    int ID;
    // we put them in hash table
    struct Hash_Node* data_item;
    struct Hash_Node* temp;
    for(int i=0; i<input_items_counter; i++)
    {
        hash_index =a[i];
        ID=hash_index;
        hash_index= con(a[i],k);
        data_item = (struct Hash_Node*)malloc(sizeof(struct Hash_Node));
        data_item->name = malloc(sizeof(char*) + 1);
        memcpy(data_item->name, names[i], sizeof(char*) + 1);   
        data_item->item = i;
        data_item->ID = ID;
        if(hash_tables[hash_index] == NULL)
        {
            hash_tables[hash_index] = data_item;
            hash_tables[hash_index]->ID = ID;
        }
        else
        {
            temp = hash_tables[hash_index];
            while(hash_tables[hash_index]->next != NULL)
            {
                hash_tables[hash_index] = hash_tables[hash_index]->next;
            }
            hash_tables[hash_index]->next = data_item;
            hash_tables[hash_index]->next->ID = ID;
            hash_tables[hash_index] = temp;
        }
    }

    int** f1 = (int**)malloc(query_items_counter * sizeof(int*));
    for(int i=0; i<query_items_counter; i++)
    {
        f1[i] = (int*)malloc(sizeof(int) * k);
    }
    for (int i=0;i<query_items_counter; i++){
        for (int j=0;j<k;j++){
            f1[i][j]=rand() % 2;
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
            h_q_result[i][j] = h_function(curves, i, dimension);   // h_function
            for(int a=0;a<i;a++){
                for(int b=0;b<j;b++){
                    if(h_q_result[i][j] == h_q_result[a][b]){   //if we find the same result we put the number 0 in f1 table//
                        f1[i][j]=0;
                    }
                }
            }
        }
    }

    int * a1=h_p(f1,k,query_items_counter);

    struct Hash_Node* temp1;
    char * nearest_neighbor =malloc(sizeof(char*)+1);
    
    for(int m=0; m<query_items_counter; m++)    // for every query show the results
    {
            hash_index= a1[m];
            int k_ID=hash_index;
            hash_index= con(a1[m],k);
            float min_dist = 1000000.0;
            temp1 = hash_tables[hash_index];
            int cp;
            hash_tables[hash_index] = temp1;
            for(int t=0;t<TableSize;t++){
                while(hash_tables[hash_index] != NULL)
                {
                    cp=0;
                    int hdi=hammingDistance(hash_index,t);
                        if((k_ID == hash_tables[hash_index]->ID || hdi<=probes ) && cp<=M)  // compare the IDs or haming distance <=probes and counter of points<=M
                    {
                    float dist = 0;
                    for(int d=1; d<=dimension; d++)
                    {
                        dist = dist + pow((query_curves[m][d] - curves[hash_tables[hash_index]->item][d]), 2);
                    }
                    dist = sqrt(dist);
                    if(dist < min_dist) // minimun LSH distance
                    {
                        min_dist = dist;
                        memcpy(nearest_neighbor, names[hash_tables[hash_index]->item], sizeof(char*) + 1);
                    }
                }
                hash_tables[hash_index] = hash_tables[hash_index]->next;        
            }
            
        }
        hash_tables[hash_index] = temp1;

        // calculate true distance and time
        float true_min_dist = 1000000.0;
        char* true_nearest_neighbor = malloc(sizeof(char*) + 1);
        for(int i=0; i<input_items_counter; i++)
        {
            float dist = 0;
            for(int d=1; d<=dimension; d++)
            {
                dist = dist + pow((query_curves[m][d] - curves[i][d]), 2);
            }
            dist = sqrt(dist);
            if(dist < true_min_dist && dist >= 0) // minimun Hypercube distance
            {
                true_min_dist = dist;
                memcpy(true_nearest_neighbor, names[i], sizeof(char*) + 1);
            }
        }

        if(min_dist != 1000000.0)
        {
            // print query
            fprintf(output_file_ptr, "Query: %s\n", query_names[m]);

            // print algorithm
            fprintf(output_file_ptr, "Algorithm: Cube_Vector\n");

            // print nearest neighbor
            fprintf(output_file_ptr, "Approximate Nearest neighbor: %s\n", nearest_neighbor);

            // print true nearest neighbor
            fprintf(output_file_ptr, "True Nearest neighbor: %s\n", true_nearest_neighbor);

            // print Hypercube distance
            fprintf(output_file_ptr, "distanceApproximate: %f\n", min_dist);

            // print true distance
            fprintf(output_file_ptr, "distanceTrue: %f\n", true_min_dist);
            // free memory
            free(true_nearest_neighbor);
        }
        else
        {
            fprintf(output_file_ptr, "Query: %s\n", query_names[m]);
            fprintf(output_file_ptr, "There is not a near vector in this bucket!\n");
        }
        fprintf(output_file_ptr, "\n");
    }


    //free memory
    for(int i=0; i<input_items_counter; i++)
    {
        free(f[i]);
        free(h_p_result[i]);
    }
    free(f);
    free(h_p_result);

    for(int i=0; i<query_items_counter; i++)
    {
        free(f1[i]);
        free(h_q_result[i]);
    }
    free(f1);
    free(h_q_result);

    free(a);
    free(a1);
    free(temp);

    free(data_item->name);
    free(data_item);
    free(temp1);
    free(nearest_neighbor);
    
}