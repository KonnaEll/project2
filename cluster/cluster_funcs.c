#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include "limits.h"
#include "cluster_funcs.h"
#include "lsh_funcs.h"
#include "cube_funcs.h"
#include "tree.h"
#include "lsh_for_frechet.h"

int min_distance(double** p, int j, double** centroid_index, int number_of_clusters, int dimension)     // minimun distance of vector from a centroid
{
    int min = INT_MAX;
    int dist;
    for(int i=0; i<number_of_clusters; i++)
    {
        dist = 0;
        for(int d=1; d<=dimension; d++)
        {
            dist = dist + pow((centroid_index[i][d] - p[j][d]), 2);
        }
        dist = sqrt(dist);
        if(dist < min)
        {
            min = dist;
        }
    }
    return min;
}

int max_distance(double* D, int input_items_counter)
{
    float max = -1;
    int max_ind;
    for(int i=0; i<input_items_counter; i++)
    {
        if(D[i] > max)
        {
            max = D[i];
            max_ind = i;
        }
    }
    return max_ind;
}

int min_distance_index(double** p, int j, double** centroid_index, int number_of_clusters, int dimension)
{
    int min = INT_MAX;
    int dist;
    int min_index;
    for(int i=0; i<number_of_clusters; i++)
    {
        dist = 0;
        for(int d=1; d<=dimension; d++)
        {
            dist = dist + pow((centroid_index[i][d] - p[j][d-1]), 2);
        }
        dist = sqrt(dist);
        // printf("%d\n", dist);
        if(dist < min)
        {
            min = dist;
            min_index = i;
            // printf("%d ", i);
        }
    }
    return min_index;
}

void classic_assign(int input_items_counter, int dimension, int number_of_clusters, double** curves, double** centroid_index, int complete, int silhouette, FILE* output_file_ptr)
{
    // Lloyd's algorithm
    int* cluster_per_item = malloc(sizeof(int) * input_items_counter);
    double** previous_centroid_index = malloc(sizeof(double*) * number_of_clusters);
    for(int i=0; i<number_of_clusters; i++)
    {
        previous_centroid_index[i] = malloc(sizeof(double) * (dimension + 1));     // save the centroids
    }
    for(int i=0; i<number_of_clusters; i++)
    {
        previous_centroid_index[i][0] = i;
        for(int z=1; z<=dimension; z++)
        {
            previous_centroid_index[i][z] = -1;
        }
    }

    struct timeval start, stop;
    gettimeofday(&start, 0);
    int count;
    while(1)
    {
        for(int j=0; j<input_items_counter; j++)
            cluster_per_item[j] = min_distance_index(curves, j, centroid_index, number_of_clusters, dimension);

        // Fix the new centroids
        int sum;
        for(int i=0; i<number_of_clusters; i++)
        {
            centroid_index[i][0] = i;
            for(int z=1; z<=dimension; z++)
            {
                count = 0;
                sum = 0;
                for(int j=0; j<input_items_counter; j++)
                {
                    if(cluster_per_item[j] == i)    // find the right centroid
                    {
                        sum = sum + curves[j][z-1];
                        count++;
                    }
                }
                if(count == 0)
                    centroid_index[i][z] = 0;
                else
                    centroid_index[i][z] = sum / count;     // average
            }
        }

        count = 0;
        for(int i=0; i<number_of_clusters; i++)     // if two centroids are the same for 2 rounds
        {
            for(int z=1; z<=dimension; z++)
            {
                if(previous_centroid_index[i][z] != centroid_index[i][z])
                {
                    count++;
                }
                previous_centroid_index[i][z] = centroid_index[i][z];
            }
        }
        if(count == 0)  // then break
        {
            break;
        }        
    }

    int j[number_of_clusters];
    for(int i=0; i<number_of_clusters; i++)
        j[i] = 0;

    int** radius = malloc(sizeof(int*) * number_of_clusters);    // array for vectors within radius
    for(int i=0; i<number_of_clusters; i++)
        radius[i] = malloc(sizeof(int) * input_items_counter);
    for(int i=0; i<number_of_clusters; i++)
        for(int j=0; j<input_items_counter; j++)
            radius[i][j] = -1;

    for(int i=0; i<number_of_clusters; i++) // fill with the items
    {
        for(int k=0; k<input_items_counter; k++)
        {
            if(cluster_per_item[k] == i)
            {
                radius[i][j[i]] = k;
                j[i]++;
            }
        }
    }

    // Clustering time
    gettimeofday(&stop, 0);
    long sec = stop.tv_sec - start.tv_sec;
    long mic_sec = stop.tv_usec - start.tv_usec;
    double clustering_time = sec + mic_sec*1e-6;


    if(complete == 0)
    {
        fprintf(output_file_ptr, "Algorithm: Lloyd's\n");   // prints
        for(int m=0; m<number_of_clusters; m++)
        {
            fprintf(output_file_ptr, "CLUSTER-%d ", m+1);
            fprintf(output_file_ptr, "(centroid: ");
            for(int i=0; i<dimension; i++)
                fprintf(output_file_ptr, "%f ", centroid_index[m][i]);
            fprintf(output_file_ptr, ")\n\n");
        }
        fprintf(output_file_ptr, "clustering time: %f\n\n", clustering_time);
    }
    else
    {
        for(int m=0; m<number_of_clusters; m++)
        {
            fprintf(output_file_ptr, "CLUSTER-%d ", m+1);
            fprintf(output_file_ptr, "{centroid: ");
            for(int i=0; i<dimension; i++)
                fprintf(output_file_ptr, "%f ", centroid_index[m][i]);
            for(int i=0; i<j[m]; i++)
                fprintf(output_file_ptr, ", %d", radius[m][i]);
            fprintf(output_file_ptr, "}\n\n");
        }
    }

    if(silhouette == 1)
    {
        // Silhouette
        long dist;
        int min;
        int min_ind[number_of_clusters];
        for(int i=0; i<number_of_clusters; i++)     // which is the closest centroid to each centroid
        {
            min = INT_MAX;
            for(int j=0; j<number_of_clusters; j++)
            {
                if(i != j)
                {
                    dist = 0;
                    for(int d=1; d<=dimension; d++)
                    {
                        dist = dist + pow((centroid_index[i][d] - centroid_index[j][d]), 2);
                    }
                    dist = sqrt(dist);
                    if(dist < min)
                    {
                        min = dist;
                        min_ind[i] = j;
                    }
                }
            }
        }

        int clust, a_i, b;
        int sum = 0;
        int s[input_items_counter];
        for(int i=0; i<input_items_counter; i++)
        {
            clust = cluster_per_item[i];
            for(int k=0; k<j[clust]-1; k++)
            {
                dist = 0;
                for(int d=1; d<=dimension; d++)
                {
                    dist = dist + pow((curves[i][d-1] - curves[radius[clust][k]][d-1]), 2);   // distance for every item
                }
                dist = sqrt(dist);

                sum = sum + dist;
            }
            a_i = sum / j[clust];   // average
            sum = 0;
            clust = min_ind[clust];     // closest centroid
            for(int k=0; k<j[clust]-1; k++)
            {
                dist = 0;
                for(int d=1; d<=dimension; d++)
                {
                    dist = dist + pow((curves[i][d-1] - curves[radius[clust][k]][d-1]), 2);
                }
                dist = sqrt(dist);
                sum = sum + dist;
            }
            b = sum / j[clust]; // average
            sum = 0;

            if(a_i > b)
                s[i] = (b / a_i) - 1;
            else if(b > a_i)
                s[i] = 1 - (a_i / b);
            else
                s[i] = 0;
        }

        float average_per_cluster[number_of_clusters];  // average for each cluster
        for(int i=0; i<number_of_clusters; i++)
            average_per_cluster[i] = 0;
        for(int i=0; i<input_items_counter; i++)
        {
            for(int m=0; m<number_of_clusters; m++)
            {
                if(cluster_per_item[i] == m)
                {
                    average_per_cluster[m] = average_per_cluster[m] + s[i];
                    continue;
                }
            }
        }

        fprintf(output_file_ptr, "Silhouette: [");
        for(int m=0; m<number_of_clusters; m++)
        {
            fprintf(output_file_ptr, "%f, ", average_per_cluster[m] / j[m]);
        }

        float total_average_sil = 0;
        for(int i=0; i<input_items_counter; i++)
            total_average_sil = total_average_sil + s[i];
        total_average_sil = total_average_sil / input_items_counter;
        fprintf(output_file_ptr, "%f]\n\n", total_average_sil);
    }


    // free memory
    free(cluster_per_item);
    for(int i=0; i<number_of_clusters; i++)
    {
        free(previous_centroid_index[i]);
        free(radius[i]);
    }
    free(previous_centroid_index);
    free(radius);
}

void lsh_assign(char** names, int input_items_counter, int dimension, int number_of_clusters, double** curves, double** centroid_index, int complete, int silhouette, FILE* output_file_ptr, int number_of_vector_hash_functions)
{
    fprintf(output_file_ptr, "Algorithm: LSH\n");
    float** h_p_result = malloc(sizeof(float*) * input_items_counter); // array with the results of the h function
    for(int i=0; i<input_items_counter; i++)
        h_p_result[i] = malloc(sizeof(float) * number_of_vector_hash_functions);
    for(int i=0; i<input_items_counter; i++)
    {
        srand(time(0));
        for(int j=0; j<number_of_vector_hash_functions; j++)
        {
            h_p_result[i][j] = h_function(curves, i, dimension);   // h_function
        }
    }

    // Hash table for input file
    int hash_index;
    int TableSize = input_items_counter / 8;
    int M = (int)pow(2, 32) - 5;
    struct Hash_Node* hash_tables[TableSize];

    for(int i=0; i<TableSize; i++)
    {
        hash_tables[i] = NULL;
    }

    // same way as lsh
    float** h_q_result = malloc(sizeof(float*) * number_of_clusters); // array with the results of the h function
    for(int i=0; i<number_of_clusters; i++)
        h_q_result[i] = malloc(sizeof(float) * number_of_vector_hash_functions);
    
    int** radius = malloc(sizeof(int*) * number_of_clusters);    // array for vectors within radius
    for(int i=0; i<number_of_clusters; i++)
        radius[i] = malloc(sizeof(int) * input_items_counter);
    
    double** previous_centroid_index = malloc(sizeof(double*) * number_of_clusters);  // save the centroids
    for(int i=0; i<number_of_clusters; i++)
        previous_centroid_index[i] = malloc(sizeof(double) * (dimension + 1));

    int r[number_of_vector_hash_functions];
    int ID;

    struct timeval start, stop;
    gettimeofday(&start, 0);
    for(int i=0; i<number_of_vector_hash_functions; i++)
    {
        r[i] = rand() % 10;
    }
    for(int i=0; i<input_items_counter; i++)
    {
        hash_index = 0;
        for(int j=0; j<number_of_vector_hash_functions; j++)
        {
            hash_index = hash_index + (int)h_p_result[i][j] * r[j];
        }
        hash_index = hash_index % M;    // mod M
        if(hash_index < 0)
        {
            hash_index = hash_index * (-1);     // only positive number
        }
        ID = hash_index;
        hash_index = hash_index % TableSize;    // mod TableSize

        struct Hash_Node* data_item = (struct Hash_Node*)malloc(sizeof(struct Hash_Node));
        data_item->item = i;
        data_item->name = malloc(sizeof(char*) + 1);
        data_item->name = names[i];
        data_item->ID = ID;
        if(hash_tables[hash_index] == NULL)
        {
            hash_tables[hash_index] = data_item;
            hash_tables[hash_index]->ID = ID;
        }
        else
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
    for(int i=0; i<number_of_clusters; i++)
    {
        for(int z=1; z<dimension; z++)
        {
            previous_centroid_index[i][z] = -1;
        }
    }

    for(int i=0; i<number_of_clusters; i++)
        for(int j=0; j<input_items_counter; j++)
            radius[i][j] = -1;

    int j[number_of_clusters];
    int j_temp[number_of_clusters];
    for(int i=0; i<number_of_clusters; i++)
    {
        j[i] = 0;
        j_temp[i] = -1;
    }

    while(1)
    {
        for(int i=0; i<number_of_clusters; i++)
        {
            srand(time(0));
            for(int j=0; j<number_of_vector_hash_functions; j++)
            {
                h_q_result[i][j] = h_function(centroid_index, i, dimension);   // h_function
            }
        }

        int min = INT_MAX;
        int dist;
        for(int i=0; i<number_of_clusters; i++)     // find first radius (minimum distance among all centroids)
        {
            for(int j=0; j<number_of_clusters; j++)
            {
                if(i != j)
                {
                    dist = 0;
                    for(int d=1; d<=dimension; d++)
                    {
                        dist = dist + pow((centroid_index[i][d] - centroid_index[j][d]), 2);
                    }
                    dist = sqrt(dist);
                    if(dist < min)
                    {
                        min = dist;
                    }
                }
            }
        }
        int R = min / 2;

        int q_hash_index[number_of_clusters];
        int k_ID[number_of_clusters];

        struct Hash_Node* temp[TableSize];
        for(int m=0; m<number_of_clusters; m++)    // for every cluster find the right bucket
        {
            q_hash_index[m] = 0;
            for(int j=0; j<number_of_vector_hash_functions; j++)
            {
                q_hash_index[m] = q_hash_index[m] + (int)h_q_result[m][j] * r[j]; // g function
            }
            q_hash_index[m] = q_hash_index[m] % M;    // mod M
            if(q_hash_index[m] < 0)
            {
                q_hash_index[m] = q_hash_index[m] * (-1);
            }
            k_ID[m] = q_hash_index[m];  // ID for comparison
            q_hash_index[m] = q_hash_index[m] % TableSize;    // mod TableSize
        }
        int count[number_of_clusters];
        for(int m=0; m<number_of_clusters; m++)
            count[m] = 0;
        int cntr = 0;

        while(1)
        {
            for(int m=0; m<number_of_clusters; m++)
            {
                temp[m] = hash_tables[q_hash_index[m]];
                while(hash_tables[q_hash_index[m]] != NULL)  // in the list
                {
                    if(k_ID[m] == hash_tables[q_hash_index[m]]->ID)  // compare the IDs
                    {
                        // calculate distance
                        int dist = 0;
                        for(int d=1; d<=dimension; d++)
                        {
                            dist = dist + pow((centroid_index[m][d] - curves[hash_tables[q_hash_index[m]]->item][d-1]), 2);
                        }
                        dist = sqrt(dist);

                        // vectors in range r
                        if(dist < R)
                        {
                            int j_dist = j[m];
                            int flag = 0;
                            for(int i=0; i<=j_dist; i++)
                            {
                                if(hash_tables[q_hash_index[m]]->item == radius[m][i])   // not same items in the list
                                {
                                    flag = 1;
                                }
                            }
                            if(flag == 0)
                            {
                                radius[m][j[m]] = hash_tables[q_hash_index[m]]->item;
                                j[m]++;
                                for(int i=0; i<number_of_clusters; i++)
                                {
                                    if(i != m)
                                    {
                                        for(int k=0; k<j[i]; k++)
                                        {
                                            // item in another cluster
                                            if(hash_tables[q_hash_index[m]]->item == radius[i][k])
                                            {
                                                int sec_dist = 0;
                                                for(int d=1; d<=dimension; d++) // calculate distance of the other
                                                {
                                                    sec_dist = sec_dist + pow((centroid_index[i][d] - curves[hash_tables[q_hash_index[m]]->item][d-1]), 2);
                                                }
                                                sec_dist = sqrt(sec_dist);
                                                if(dist < sec_dist) // closer cluster
                                                {
                                                    radius[i][k] = -1;
                                                }
                                                else
                                                {
                                                    radius[m][j[m] - 1] = -1;
                                                    j[m]--;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    hash_tables[q_hash_index[m]] = hash_tables[q_hash_index[m]]->next;
                }
                hash_tables[q_hash_index[m]] = temp[m];
                if(j_temp[m] == j[m])   // if no item gets in any cluster
                {
                    count[m]++;
                }
                j_temp[m] = j[m];
            }
            for(int i=0; i<number_of_clusters; i++)
            {
                if(count[i] >= 3)   // for 3 times
                {
                    cntr++;
                }
            }
            if(cntr > 0)    // then break
                break;
            R = R * 2;  // or double the radius and continue
        }

        // Fix the new centroids
        int sum;
        for(int i=0; i<number_of_clusters; i++)
        {
            centroid_index[i][0] = i;
            for(int z=1; z<dimension; z++)
            {
                sum = 0;
                for(int k=0; k<j[i]; k++)
                {
                    sum = sum + curves[radius[i][k]][z-1]; // average
                }
                if(j[i] > 0)
                    centroid_index[i][z] = sum / j[i];
            }
        }

        int counter = 0;
        for(int i=0; i<number_of_clusters; i++)     // centroids are the same for two rounds
        {
            for(int z=1; z<dimension; z++)
            {
                if(previous_centroid_index[i][z] != centroid_index[i][z])
                {
                    counter++;
                }
                previous_centroid_index[i][z] = centroid_index[i][z];
            }
        }

        if(counter == 0)    // then stop
        {
            break;
        }
    }

    int count = 0;
    int min = -1;
    int centroids[input_items_counter]; // cluster of each item
    for(int k=0; k<input_items_counter; k++)
    {
        for(int m=0; m<number_of_clusters; m++)
        {
            for(int i=0; i<j[m]; i++)
            {
                if(curves[k][0] == radius[m][i])
                {
                    centroids[k] = m;
                    count++;
                }
            }
        }
        if(count == 0)
        {
            // whoever is left behind then find the closer cluster and put it in there
            min = min_distance_index(curves, k, centroid_index, number_of_clusters, dimension);
            radius[min][j[min]] = k;
            j[min]++;
            centroids[k] = min;
        }
        count = 0;
    }

    // Clustering time
    gettimeofday(&stop, 0);
    long sec = stop.tv_sec - start.tv_sec;
    long mic_sec = stop.tv_usec - start.tv_usec;
    double clustering_time = sec + mic_sec*1e-6;

    if(complete == 1)
    {
        for(int m=0; m<number_of_clusters; m++)
        {
            fprintf(output_file_ptr, "CLUSTER-%d ", m+1);
            fprintf(output_file_ptr, "{centroid: ");
            for(int i=0; i<dimension; i++)
                fprintf(output_file_ptr, "%f ", centroid_index[m][i]);
            for(int i=0; i<j[m]; i++)
                fprintf(output_file_ptr, ", %d", radius[m][i]);
            fprintf(output_file_ptr, "}\n\n");
        }
    }
    else
    {
        for(int m=0; m<number_of_clusters; m++)
        {
            fprintf(output_file_ptr, "CLUSTER-%d ", m+1);
            fprintf(output_file_ptr, "{size: %d, ", j[m]);
            fprintf(output_file_ptr, "centroid: ");
            for(int i=0; i<dimension; i++)
                fprintf(output_file_ptr, "%f ", centroid_index[m][i]);
            fprintf(output_file_ptr, "}\n\n");
        }
        fprintf(output_file_ptr, "clustering time: %f\n\n", clustering_time);
    }
    if(silhouette == 1)
    {
        // Silhouette
        long dist;
        int min_ind[number_of_clusters];
        for(int i=0; i<number_of_clusters; i++)     // which is the closest centroid to each centroid
        {
            min = INT_MAX;
            for(int j=0; j<number_of_clusters; j++)
            {
                if(i != j)
                {
                    dist = 0;
                    for(int d=1; d<=dimension; d++)
                    {
                        dist = dist + pow((centroid_index[i][d] - centroid_index[j][d]), 2);
                    }
                    dist = sqrt(dist);
                    if(dist < min)
                    {
                        min = dist;
                        min_ind[i] = j;
                    }
                }
            }
        }

        int clust, a_i, b;
        double sum = 0;
        int s[input_items_counter];
        for(int i=0; i<input_items_counter; i++)
        {
            clust = centroids[i];
            for(int k=0; k<j[clust]-1; k++)
            {
                dist = 0;
                for(int d=1; d<=dimension; d++)
                {
                    dist = dist + pow((curves[i][d-1] - curves[radius[clust][k]][d-1]), 2);   // distance for every item
                }
                dist = sqrt(dist);
                sum = sum + dist;
            }
            if(j[clust] == 0)
                a_i = 1;
            else
                a_i = sum / j[clust];   // average
            
            sum = 0;
            clust = min_ind[clust];     // closest centroid
            for(int k=0; k<j[clust]-1; k++)
            {
                dist = 0;
                for(int d=1; d<=dimension; d++)
                {
                    dist = dist + pow((curves[i][d-1] - curves[radius[clust][k]][d-1]), 2);
                }
                dist = sqrt(dist);
                sum = sum + dist;
            }
            if(j[clust] == 0)
                b = 1;
            else
                b = sum / j[clust]; // average
            sum = 0;

            if(a_i > b)
                s[i] = (b / a_i) - 1;
            else if(b > a_i)
                s[i] = 1 - (a_i / b);
            else
                s[i] = 0;
        }

        float average_per_cluster[number_of_clusters];  // average for each cluster
        for(int i=0; i<number_of_clusters; i++)
            average_per_cluster[i] = 0;
        for(int i=0; i<input_items_counter; i++)
        {
            for(int m=0; m<number_of_clusters; m++)
            {
                if(centroids[i] == m)
                {
                    average_per_cluster[m] = average_per_cluster[m] + s[i];
                    continue;
                }
            }
        }

        fprintf(output_file_ptr, "Silhouette: [");
        for(int m=0; m<number_of_clusters; m++)
        {
            fprintf(output_file_ptr, "%f, ", average_per_cluster[m] / j[m]);
        }

        float total_average_sil = 0;
        for(int i=0; i<input_items_counter; i++)
            total_average_sil = total_average_sil + s[i];
        total_average_sil = total_average_sil / input_items_counter;
        fprintf(output_file_ptr, "%f]\n\n\n", total_average_sil);
    }

    // free memory
    for(int i=0; i<input_items_counter; i++)
    {
        free(h_p_result[i]);
    }
    free(h_p_result);
    for(int i=0; i<number_of_clusters; i++)
    {
        free(previous_centroid_index[i]);
        // free(radius[i]);
        free(h_q_result[i]);
    }
    free(previous_centroid_index);
    // free(radius);
    free(h_q_result);
}


void cube_assign(char** names,int input_items_counter, int dimension, int number_of_clusters, double** curves, double** centroid_index, int complete, int silhouette, FILE* output_file_ptr,int k, int probes, int M)
{
    // create 2 tables with 2 dimension and we put inside 0,1//
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

    int** f1 = (int**)malloc(number_of_clusters * sizeof(int*));
    for(int i=0; i<number_of_clusters; i++)
    {
        f1[i] = (int*)malloc(sizeof(int) * k);
    }
    for (int i=0;i<number_of_clusters; i++){
        for (int  j=0;j<k;j++){
            f1[i][j]=rand() % 2;
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

    // Hash table for input file same as before
    int TableSize = pow(2,k);
    int hash_index;
    struct Hash_Node* hash_tables[TableSize];

    for(int i=0; i<(TableSize); i++)
    {
        hash_tables[i] = NULL;
    }

    int ID;
    struct timeval start, stop;
    gettimeofday(&start, 0);
    for(int i=0; i<input_items_counter; i++)
    {
        hash_index =a[i];
        ID=hash_index;
        hash_index= con(a[i],k);
        struct Hash_Node* data_item = (struct Hash_Node*)malloc(sizeof(struct Hash_Node));
        data_item->item = i;
        data_item->name = malloc(sizeof(char*) + 1);
        data_item->name = names[i];
        data_item->ID = ID;
        if(hash_tables[hash_index] == NULL)
        {
            hash_tables[hash_index] = data_item;
            hash_tables[hash_index]->ID = ID;
        }
        else
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
    int** previous_centroid_index = malloc(sizeof(int*) * number_of_clusters);
    for(int i=0; i<number_of_clusters; i++)
    {
        previous_centroid_index[i] = malloc(sizeof(int) * (dimension + 1));
    }
    for(int i=0; i<number_of_clusters; i++)
    {
        for(int z=1; z<dimension; z++)
        {
            previous_centroid_index[i][z] = -1;
        }
    }
    int** radius = malloc(sizeof(int*) * number_of_clusters);    // array for vectors within radius
    for(int i=0; i<number_of_clusters; i++)
        radius[i] = malloc(sizeof(int) * input_items_counter);
    for(int i=0; i<number_of_clusters; i++)
        for(int j=0; j<input_items_counter; j++)
            radius[i][j] = -1;

    int j[number_of_clusters];
    int j_temp[number_of_clusters];
    for(int i=0; i<number_of_clusters; i++)
    {
        j[i] = 0;
        j_temp[i] = -1;
    }

    float** h_q_result = malloc(sizeof(float*) * number_of_clusters); // array with the results of the h function
    for(int i=0; i<number_of_clusters; i++)
        h_q_result[i] = malloc(sizeof(float) * k);

    while(1)
    {
        for(int i=0; i<number_of_clusters; i++)
        {
            srand(time(0));
            for(int j=0; j<k; j++)
            {
                h_q_result[i][j] = h_function(centroid_index, i, dimension);   // h_function
                for(int a=0;a<i;a++){
                    for(int b=0;b<j;b++){
                        if(h_q_result[i][j] == h_q_result[a][b]){   //if we find the same result we put the number 0 in f1 table//
                            f1[i][j]=0;
                        }
                    }
                }
            }
        }

        int min = INT_MAX;
        long dist;                        
        for(int i=0; i<number_of_clusters; i++)     // find first radius (minimum distance among all centroids)
        {
            for(int j=0; j<number_of_clusters; j++)
            {
                if(i != j)
                {
                    dist = 0;
                    for(int d=1; d<=dimension; d++)
                    {
                        dist = dist + pow((centroid_index[i][d] - centroid_index[j][d]), 2);
                    }
                    dist = sqrt(dist);
                    if(dist < min)
                    {
                        min = dist;
                    }
                }
            }
        }

        int R = min / 2;
        
        struct Hash_Node* temp[TableSize];
        int q_hash_index[number_of_clusters];
        int k_ID[number_of_clusters];
        int* a1=h_p(f1,k,number_of_clusters);
        for(int m=0; m<number_of_clusters; m++)    // for every query show the results
        {
            q_hash_index[m]=a1[m];
            k_ID[m]=q_hash_index[m];
            q_hash_index[m]=con(a1[m],k);
        }                
        int count[number_of_clusters];
        for(int m=0; m<number_of_clusters; m++)
            count[m] = 0;
        int cntr = 0;
        int cp;
        // exactly the same as lsh before but with the cube algorithm
        while(1)
        {
            for(int m=0; m<number_of_clusters; m++)    // for every query show the results
            {
                temp[m] = hash_tables[q_hash_index[m]];
                for(int t=0;t<TableSize;t++)
                {
                    while(hash_tables[q_hash_index[m]] != NULL)
                    {
                        cp=0;
                        // int hdi=hammingDistance(hash_index,t);
                        // || hdi<=probes

                        if((k_ID[m] == hash_tables[q_hash_index[m]]->ID ) && cp<=M)  // compare the IDs
                        {
                            // calculate distance
                            dist = 0;
                            for(int d=1; d<=dimension; d++)
                            {
                                dist = dist + pow((centroid_index[m][d] - curves[hash_tables[q_hash_index[m]]->item][d-1]), 2);
                            }
                            dist = sqrt(dist);
                            // vectors in range r
                            if(dist < R)
                            {
                                int j_dist = j[m];
                                int flag = 0;
                                for(int i=0; i<=j_dist; i++)
                                {
                                    if(hash_tables[q_hash_index[m]]->item == radius[m][i])
                                    {
                                        flag = 1;
                                    }
                                }
                                if(flag == 0)
                                {
                                    radius[m][j[m]] = hash_tables[q_hash_index[m]]->item;
                                    j[m]++;
                                    for(int i=0; i<number_of_clusters; i++)
                                    {
                                        if(i != m)
                                        {
                                            for(int k=0; k<j[i]; k++)
                                            {
                                                if(hash_tables[q_hash_index[m]]->item == radius[i][k])
                                                {
                                                    long sec_dist = 0;
                                                    for(int d=1; d<=dimension; d++)
                                                    {
                                                        sec_dist = sec_dist + pow((centroid_index[i][d] - curves[hash_tables[q_hash_index[m]]->item][d-1]), 2);
                                                    }
                                                    sec_dist = sqrt(sec_dist);
                                                    if(dist < sec_dist)
                                                    {
                                                        radius[i][k] = -1;
                                                    }
                                                    else
                                                    {
                                                        radius[m][j[m] - 1] = -1;
                                                        j[m]--;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        hash_tables[q_hash_index[m]] = hash_tables[q_hash_index[m]]->next;
                    }
                }
                hash_tables[q_hash_index[m]] = temp[m];
                if(j_temp[m] == j[m])   // when to break
                {
                    count[m]++;
                }
                j_temp[m] = j[m];
            }
            for(int i=0; i<number_of_clusters; i++)
            {
                if(count[i] >= 3)
                {
                    cntr++;
                }
            }
            if(cntr > 0)
                break;
            R = R * 2;
        }

        // Fix the new centroids
        int sum;
        for(int i=0; i<number_of_clusters; i++)
        {
            centroid_index[i][0] = i;
            for(int z=1; z<dimension; z++)
            {
                sum = 0;
                for(int k=0; k<j[i]; k++)
                {
                    sum = sum + curves[radius[i][k]][z-1];
                }
                if(j[i] > 0)
                    centroid_index[i][z] = sum / j[i];
            }
        }

        int counter = 0;
        for(int i=0; i<number_of_clusters; i++)
        {
            for(int z=1; z<dimension; z++)
            {
                if(previous_centroid_index[i][z] != centroid_index[i][z])
                {
                    counter++;
                }
                previous_centroid_index[i][z] = centroid_index[i][z];
            }
        }

        if(counter == 0)
        {
            break;
        }
    }
        
    int count = 0;
    int min = -1;
    int centroids[input_items_counter];
    for(int k=0; k<input_items_counter; k++)
    {
        for(int m=0; m<number_of_clusters; m++)
        {
            for(int i=0; i<j[m]; i++)
            {
                if(curves[k][0] == radius[m][i])
                {
                    centroids[k] = m;
                    count++;
                }
            }
        }
        if(count == 0)
        {
            min = min_distance_index(curves, k, centroid_index, number_of_clusters, dimension);  // what is left
            radius[min][j[min]] = k;
            j[min]++;
            centroids[k] = min;
        }
        count = 0;
    }
    // Clustering time
    gettimeofday(&stop, 0);
    long sec = stop.tv_sec - start.tv_sec;
    long mic_sec = stop.tv_usec - start.tv_usec;
    double clustering_time = sec + mic_sec*1e-6;

    if(complete == 1)
    {   for(int m=0; m<number_of_clusters; m++)
        {
            fprintf(output_file_ptr, "CLUSTER-%d ", m+1);
            fprintf(output_file_ptr, "centroid: ");
            for(int i=0; i<dimension; i++)
                fprintf(output_file_ptr, "%f ", centroid_index[m][i]);
            for(int i=0; i<j[m]; i++)
                fprintf(output_file_ptr, ", %d", radius[m][i]);
            fprintf(output_file_ptr, ")\n\n");
        }
    }
    else
    {
        for(int m=0; m<number_of_clusters; m++)
        {
            fprintf(output_file_ptr, "CLUSTER-%d ", m+1);
            fprintf(output_file_ptr, "size: %d, ", j[m]);
            fprintf(output_file_ptr, "(centroid: ");
            for(int i=0; i<dimension; i++)
                fprintf(output_file_ptr, "%f ", centroid_index[m][i]);
            fprintf(output_file_ptr, ")\n\n");
        }
        fprintf(output_file_ptr, "clustering time: %f\n\n", clustering_time);
    }
    
    if(silhouette == 1)
    {
        // Silhouette is exactly the same
        long dist;
        int min_ind[number_of_clusters];
        for(int i=0; i<number_of_clusters; i++)     // which is the closest centroid to each centroid
        {
            min = INT_MAX;
            for(int j=0; j<number_of_clusters; j++)
            {
                if(i != j)
                {
                    dist = 0;
                    for(int d=1; d<=dimension; d++)
                    {
                        dist = dist + pow((centroid_index[i][d] - centroid_index[j][d]), 2);
                    }
                    dist = sqrt(dist);
                    if(dist < min)
                    {
                        min = dist;
                        min_ind[i] = j;
                    }
                }
            }
        }

        int clust, a_i, b;
        double sum = 0;
        int s[input_items_counter];
        for(int i=0; i<input_items_counter; i++)
        {
            clust = centroids[i];
            for(int k=0; k<j[clust]-1; k++)
            {
                dist = 0;
                for(int d=1; d<=dimension; d++)
                {
                    dist = dist + pow((curves[i][d-1] - curves[radius[clust][k]][d-1]), 2);
                }
                dist = sqrt(dist);
                sum = sum + dist;
            }
            if(j[clust] == 0)
                a_i = 1;
            else
                a_i = sum / j[clust];   // average;
            sum = 0;

            clust = min_ind[clust];     // closest centroid
            for(int k=0; k<j[clust]-1; k++)
            {
                dist = 0;
                for(int d=1; d<=dimension; d++)
                {
                    dist = dist + pow((curves[i][d-1] - curves[radius[clust][k]][d-1]), 2);
                }
                dist = sqrt(dist);
                sum = sum + dist;
            }
            if(j[clust] == 0)
                b = 1;
            else
                b = sum / j[clust]; // average
            sum = 0;

            if(a_i > b)
                s[i] = (b / a_i) -1;
            else if(b > a_i)
                s[i] = 1 - (a_i / b);
            else
                s[i] = 0;
        }

        float average_per_cluster[number_of_clusters];
        for(int i=0; i<number_of_clusters; i++)
            average_per_cluster[i] = 0;
        for(int i=0; i<input_items_counter; i++)
        {
            for(int m=0; m<number_of_clusters; m++)
            {
                if(centroids[i] == m)
                {
                    average_per_cluster[m] = average_per_cluster[m] + s[i];
                    continue;
                }
            }
        }

        fprintf(output_file_ptr, "Silhouette: [");
        for(int m=0; m<number_of_clusters; m++)
        {
            if(j[m] == 0)
                printf("00000\n");
            fprintf(output_file_ptr, "%f, ", average_per_cluster[m] / j[m]);
        }

        float total_average_sil = 0;
        for(int i=0; i<input_items_counter; i++)
            total_average_sil = total_average_sil + s[i];
        total_average_sil = total_average_sil / input_items_counter;
        fprintf(output_file_ptr, "%f]\n\n", total_average_sil);
    }

    // free memory
    for(int i=0; i<input_items_counter; i++)
    {
        free(h_p_result[i]);
    }
    free(h_p_result);
    for(int i=0; i<number_of_clusters; i++)
    {
        free(previous_centroid_index[i]);
        // free(radius[i]);
        free(h_q_result[i]);
    }
    free(previous_centroid_index);
    // free(radius);
    free(h_q_result);
}

void frechet_classic_assign(int input_items_counter, int dimension, int number_of_clusters, double** curves, double** centroid_index, int complete, int silhouette, FILE* output_file_ptr)
{
    // Lloyd's algorithm
    int* cluster_per_item = malloc(sizeof(int) * input_items_counter);
    double** previous_centroid_index = malloc(sizeof(double*) * number_of_clusters);
    for(int i=0; i<number_of_clusters; i++)
    {
        previous_centroid_index[i] = malloc(sizeof(double) * (dimension + 1));     // save the centroids
    }
    for(int i=0; i<number_of_clusters; i++)
    {
        previous_centroid_index[i][0] = i;
        for(int z=1; z<=dimension; z++)
        {
            previous_centroid_index[i][z] = -1;
        }
    }

    double* centre = malloc(sizeof(double*) * dimension);


    // int h = ceil(log2(input_items_counter));
    // int height = pow(2, h);

    struct timeval start, stop;
    gettimeofday(&start, 0);
    int count;
    while(1)
    {
        printf("ee\n");
        for(int j=0; j<input_items_counter; j++)
            cluster_per_item[j] = min_distance_index(curves, j, centroid_index, number_of_clusters, dimension);

        // Fix the new centroids
        srand((unsigned int)time(NULL));
        for(int i=0; i<number_of_clusters; i++)
        {
            centroid_index[i][0] = i;

            struct Tree_Node* root = NULL;
            struct Queue* queue = New_Queue(input_items_counter);
            for(int i=0; i<input_items_counter; i++)
            {
                Insert_node(&root, curves[rand() % input_items_counter], queue);
            }

            centre = post_order_traversal(root, dimension);
            for(int d=1; d<dimension; d++)
            {
                centroid_index[i][d] = centre[d];
            }
        }

        count = 0;
        for(int i=0; i<number_of_clusters; i++)     // if two centroids are the same for 2 rounds
        {
            for(int z=1; z<=dimension; z++)
            {
                if(previous_centroid_index[i][z] != centroid_index[i][z])
                {
                    count++;
                }
                previous_centroid_index[i][z] = centroid_index[i][z];
            }
        }
        if(count == 0)  // then break
        {
            break;
        }
    }
    // for(int i=0; i<number_of_clusters; i++)
    // {
    //     printf("\n\n");
    //     for(int j=0; j<dimension; j++)
    //         printf("%f ", centroid_index[i][j]);
    // }

    int j[number_of_clusters];
    for(int i=0; i<number_of_clusters; i++)
        j[i] = 0;

    int** radius = malloc(sizeof(int*) * number_of_clusters);    // array for vectors within radius
    for(int i=0; i<number_of_clusters; i++)
        radius[i] = malloc(sizeof(int) * input_items_counter);
    for(int i=0; i<number_of_clusters; i++)
        for(int j=0; j<input_items_counter; j++)
            radius[i][j] = -1;

    for(int i=0; i<number_of_clusters; i++) // fill with the items
    {
        for(int k=0; k<input_items_counter; k++)
        {
            if(cluster_per_item[k] == i)
            {
                radius[i][j[i]] = k;
                j[i]++;
            }
        }
    }

    // Clustering time
    gettimeofday(&stop, 0);
    long sec = stop.tv_sec - start.tv_sec;
    long mic_sec = stop.tv_usec - start.tv_usec;
    double clustering_time = sec + mic_sec*1e-6;


    if(complete == 0)
    {
        fprintf(output_file_ptr, "Algorithm: Lloyd's\n");   // prints
        for(int m=0; m<number_of_clusters; m++)
        {
            fprintf(output_file_ptr, "CLUSTER-%d ", m+1);
            fprintf(output_file_ptr, "(centroid: ");
            for(int i=0; i<dimension; i++)
                fprintf(output_file_ptr, "%f ", centroid_index[m][i]);
            fprintf(output_file_ptr, ")\n\n");
        }
        fprintf(output_file_ptr, "clustering time: %f\n\n", clustering_time);
    }
    else
    {
        for(int m=0; m<number_of_clusters; m++)
        {
            fprintf(output_file_ptr, "CLUSTER-%d ", m+1);
            fprintf(output_file_ptr, "{centroid: ");
            for(int i=0; i<dimension; i++)
                fprintf(output_file_ptr, "%f ", centroid_index[m][i]);
            for(int i=0; i<j[m]; i++)
                fprintf(output_file_ptr, ", %d", radius[m][i]);
            fprintf(output_file_ptr, "}\n\n");
        }
    }

    if(silhouette == 1)
    {
        // Silhouette
        long dist;
        int min;
        int min_ind[number_of_clusters];
        for(int i=0; i<number_of_clusters; i++)     // which is the closest centroid to each centroid
        {
            min = INT_MAX;
            for(int j=0; j<number_of_clusters; j++)
            {
                if(i != j)
                {
                    dist = 0;
                    for(int d=1; d<=dimension; d++)
                    {
                        dist = dist + pow((centroid_index[i][d] - centroid_index[j][d]), 2);
                    }
                    dist = sqrt(dist);
                    if(dist < min)
                    {
                        min = dist;
                        min_ind[i] = j;
                    }
                }
            }
        }

        int clust, a_i, b;
        int sum = 0;
        int s[input_items_counter];
        for(int i=0; i<input_items_counter; i++)
        {
            clust = cluster_per_item[i];
            for(int k=0; k<j[clust]-1; k++)
            {
                dist = 0;
                for(int d=1; d<=dimension; d++)
                {
                    dist = dist + pow((curves[i][d-1] - curves[radius[clust][k]][d-1]), 2);   // distance for every item
                }
                dist = sqrt(dist);

                sum = sum + dist;
            }
            a_i = sum / j[clust];   // average
            sum = 0;
            clust = min_ind[clust];     // closest centroid
            for(int k=0; k<j[clust]-1; k++)
            {
                dist = 0;
                for(int d=1; d<=dimension; d++)
                {
                    dist = dist + pow((curves[i][d-1] - curves[radius[clust][k]][d-1]), 2);
                }
                dist = sqrt(dist);
                sum = sum + dist;
            }
            b = sum / j[clust]; // average
            sum = 0;

            if(a_i > b)
                s[i] = (b / a_i) - 1;
            else if(b > a_i)
                s[i] = 1 - (a_i / b);
            else
                s[i] = 0;
        }

        float average_per_cluster[number_of_clusters];  // average for each cluster
        for(int i=0; i<number_of_clusters; i++)
            average_per_cluster[i] = 0;
        for(int i=0; i<input_items_counter; i++)
        {
            for(int m=0; m<number_of_clusters; m++)
            {
                if(cluster_per_item[i] == m)
                {
                    average_per_cluster[m] = average_per_cluster[m] + s[i];
                    continue;
                }
            }
        }

        fprintf(output_file_ptr, "Silhouette: [");
        for(int m=0; m<number_of_clusters; m++)
        {
            fprintf(output_file_ptr, "%f, ", average_per_cluster[m] / j[m]);
        }

        float total_average_sil = 0;
        for(int i=0; i<input_items_counter; i++)
            total_average_sil = total_average_sil + s[i];
        total_average_sil = total_average_sil / input_items_counter;
        fprintf(output_file_ptr, "%f]\n\n", total_average_sil);
    }


    // free memory
    free(cluster_per_item);
    for(int i=0; i<number_of_clusters; i++)
    {
        free(previous_centroid_index[i]);
        free(radius[i]);
    }
    free(previous_centroid_index);
    free(radius);
}


void frechet_lsh_assign(char** names, int input_items_counter, int dimension, int number_of_clusters, double** curves, double** centroid_index, int complete, int silhouette, FILE* output_file_ptr, int number_of_vector_hash_functions)
{
    fprintf(output_file_ptr, "Algorithm: LSH with Frechet\n");
    float** h_p_result = malloc(sizeof(float*) * input_items_counter); // array with the results of the h function
    for(int i=0; i<input_items_counter; i++)
        h_p_result[i] = malloc(sizeof(float) * number_of_vector_hash_functions);
    for(int i=0; i<input_items_counter; i++)
    {
        srand(time(0));
        for(int j=0; j<number_of_vector_hash_functions; j++)
        {
            h_p_result[i][j] = h_function(curves, i, dimension);   // h_function
        }
    }

    // Hash table for input file
    int hash_index;
    int TableSize = input_items_counter / 8;
    int M = (int)pow(2, 32) - 5;
    struct Hash_Node* hash_tables[TableSize];

    for(int i=0; i<TableSize; i++)
    {
        hash_tables[i] = NULL;
    }

    // same way as lsh
    float** h_q_result = malloc(sizeof(float*) * number_of_clusters); // array with the results of the h function
    for(int i=0; i<number_of_clusters; i++)
        h_q_result[i] = malloc(sizeof(float) * number_of_vector_hash_functions);
    
    int** radius = malloc(sizeof(int*) * number_of_clusters);    // array for vectors within radius
    for(int i=0; i<number_of_clusters; i++)
        radius[i] = malloc(sizeof(int) * input_items_counter);
    
    double** previous_centroid_index = malloc(sizeof(double*) * number_of_clusters);  // save the centroids
    for(int i=0; i<number_of_clusters; i++)
        previous_centroid_index[i] = malloc(sizeof(double) * (dimension + 1));

    int r[number_of_vector_hash_functions];
    int ID;

    struct timeval start, stop;
    gettimeofday(&start, 0);
    for(int i=0; i<number_of_vector_hash_functions; i++)
    {
        r[i] = rand() % 10;
    }
    for(int i=0; i<input_items_counter; i++)
    {
        hash_index = 0;
        for(int j=0; j<number_of_vector_hash_functions; j++)
        {
            hash_index = hash_index + (int)h_p_result[i][j] * r[j];
        }
        hash_index = hash_index % M;    // mod M
        if(hash_index < 0)
        {
            hash_index = hash_index * (-1);     // only positive number
        }
        ID = hash_index;
        hash_index = hash_index % TableSize;    // mod TableSize

        struct Hash_Node* data_item = (struct Hash_Node*)malloc(sizeof(struct Hash_Node));
        data_item->item = i;
        data_item->name = malloc(sizeof(char*) + 1);
        data_item->name = names[i];
        data_item->ID = ID;
        if(hash_tables[hash_index] == NULL)
        {
            hash_tables[hash_index] = data_item;
            hash_tables[hash_index]->ID = ID;
        }
        else
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
    for(int i=0; i<number_of_clusters; i++)
    {
        for(int z=1; z<dimension; z++)
        {
            previous_centroid_index[i][z] = -1;
        }
    }

    for(int i=0; i<number_of_clusters; i++)
        for(int j=0; j<input_items_counter; j++)
            radius[i][j] = -1;

    int j[number_of_clusters];
    int j_temp[number_of_clusters];
    for(int i=0; i<number_of_clusters; i++)
    {
        j[i] = 0;
        j_temp[i] = -1;
    }

    double* centre = malloc(sizeof(double*) * dimension);
    int counter = 0;
    float** distance = malloc(sizeof(float*) * dimension);
    for(int i=0; i<dimension; i++)
        distance[i] = malloc(sizeof(float) * dimension);
    while(1)
    {
        counter++;
        for(int i=0; i<number_of_clusters; i++)
        {
            srand(time(0));
            for(int j=0; j<number_of_vector_hash_functions; j++)
            {
                h_q_result[i][j] = h_function(centroid_index, i, dimension);   // h_function
            }
        }

        int min = INT_MAX;
        int dist;
        for(int i=0; i<number_of_clusters; i++)     // find first radius (minimum distance among all centroids)
        {
            for(int j=0; j<number_of_clusters; j++)
            {
                if(i != j)
                {
                    dist = 0;
                    for(int d=1; d<=dimension; d++)
                    {
                        dist = dist + pow((centroid_index[i][d] - centroid_index[j][d]), 2);
                    }
                    dist = sqrt(dist);
                    if(dist < min)
                    {
                        min = dist;
                    }
                }
            }
        }
        int R = min / 2;

        int q_hash_index[number_of_clusters];
        int k_ID[number_of_clusters];

        struct Hash_Node* temp[TableSize];
        for(int m=0; m<number_of_clusters; m++)    // for every cluster find the right bucket
        {
            q_hash_index[m] = 0;
            for(int j=0; j<number_of_vector_hash_functions; j++)
            {
                q_hash_index[m] = q_hash_index[m] + (int)h_q_result[m][j] * r[j]; // g function
            }
            q_hash_index[m] = q_hash_index[m] % M;    // mod M
            if(q_hash_index[m] < 0)
            {
                q_hash_index[m] = q_hash_index[m] * (-1);
            }
            k_ID[m] = q_hash_index[m];  // ID for comparison
            q_hash_index[m] = q_hash_index[m] % TableSize;    // mod TableSize
        }
        int count[number_of_clusters];
        for(int m=0; m<number_of_clusters; m++)
            count[m] = 0;
        int cntr = 0;

        while(1)
        {
            for(int m=0; m<number_of_clusters; m++)
            {
                temp[m] = hash_tables[q_hash_index[m]];
                while(hash_tables[q_hash_index[m]] != NULL)  // in the list
                {
                    if(k_ID[m] == hash_tables[q_hash_index[m]]->ID)  // compare the IDs
                    {
                        // // calculate distance
                        // int dist = 0;
                        // for(int d=1; d<=dimension; d++)
                        // {
                        //     dist = dist + pow((centroid_index[m][d] - curves[hash_tables[q_hash_index[m]]->item][d-1]), 2);
                        // }
                        // dist = sqrt(dist);

                        // calculate distance
                        for(int i=0; i<dimension; i++)
                            for(int j=0; j<dimension; j++)
                                distance[i][j] = -1;
                        dist = distance_computation(distance, dimension, curves, centroid_index, m, hash_tables[q_hash_index[m]]->item);

                        // vectors in range r
                        if(dist < R)
                        {
                            int j_dist = j[m];
                            int flag = 0;
                            for(int i=0; i<=j_dist; i++)
                            {
                                if(hash_tables[q_hash_index[m]]->item == radius[m][i])   // not same items in the list
                                {
                                    flag = 1;
                                }
                            }
                            if(flag == 0)
                            {
                                radius[m][j[m]] = hash_tables[q_hash_index[m]]->item;
                                j[m]++;
                                for(int i=0; i<number_of_clusters; i++)
                                {
                                    if(i != m)
                                    {
                                        for(int k=0; k<j[i]; k++)
                                        {
                                            // item in another cluster
                                            if(hash_tables[q_hash_index[m]]->item == radius[i][k])
                                            {
                                                int sec_dist = 0;
                                                // for(int d=1; d<=dimension; d++) // calculate distance of the other
                                                // {
                                                //     sec_dist = sec_dist + pow((centroid_index[i][d] - curves[hash_tables[q_hash_index[m]]->item][d-1]), 2);
                                                // }
                                                // sec_dist = sqrt(sec_dist);
                                                for(int i=0; i<dimension; i++)
                                                    for(int j=0; j<dimension; j++)
                                                        distance[i][j] = -1;
                                                sec_dist = distance_computation(distance, dimension, curves, centroid_index, i, hash_tables[q_hash_index[m]]->item);

                                                if(dist < sec_dist) // closer cluster
                                                {
                                                    radius[i][k] = -1;
                                                }
                                                else
                                                {
                                                    radius[m][j[m] - 1] = -1;
                                                    j[m]--;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    hash_tables[q_hash_index[m]] = hash_tables[q_hash_index[m]]->next;
                }
                hash_tables[q_hash_index[m]] = temp[m];
                if(j_temp[m] == j[m])   // if no item gets in any cluster
                {
                    count[m]++;
                }
                j_temp[m] = j[m];
            }
            for(int i=0; i<number_of_clusters; i++)
            {
                if(count[i] >= 3)   // for 3 times
                {
                    cntr++;
                }
            }
            if(cntr > 0)    // then break
                break;
            R = R * 2;  // or double the radius and continue
        }
        // Fix the new centroids
        srand((unsigned int)time(NULL));
        for(int i=0; i<number_of_clusters; i++)
        {
            centroid_index[i][0] = i;

            struct Tree_Node* root = NULL;
            struct Queue* queue = New_Queue(input_items_counter);
            for(int i=0; i<input_items_counter; i++)
            {
                Insert_node(&root, curves[rand() % input_items_counter], queue);
            }

            centre = post_order_traversal(root, dimension);
            for(int d=1; d<dimension; d++)
            {
                centroid_index[i][d] = centre[d];
            }
        }

        // counter = 0;
        // for(int i=0; i<number_of_clusters; i++)     // if two centroids are the same for 2 rounds
        // {
        //     for(int z=1; z<dimension; z++)
        //     {
        //         if(previous_centroid_index[i][z] != centroid_index[i][z])
        //         {
        //             counter++;
        //         }
        //         previous_centroid_index[i][z] = centroid_index[i][z];
        //     }
        // }
        if(counter > 10)  // then break
        {
            break;
        }
        // for(int i=0; i<number_of_clusters; i++)
        // {
        //     printf("\n\n");
        //     for(int j=0; j<dimension; j++)
        //         printf("%f ", centroid_index[i][j]);
        // }
        // exit(1);
    }

    int count = 0;
    int min = -1;
    int centroids[input_items_counter]; // cluster of each item
    for(int k=0; k<input_items_counter; k++)
    {
        for(int m=0; m<number_of_clusters; m++)
        {
            for(int i=0; i<j[m]; i++)
            {
                if(curves[k][0] == radius[m][i])
                {
                    centroids[k] = m;
                    count++;
                }
            }
        }
        if(count == 0)
        {
            // whoever is left behind then find the closer cluster and put it in there
            min = min_distance_index(curves, k, centroid_index, number_of_clusters, dimension);
            radius[min][j[min]] = k;
            j[min]++;
            centroids[k] = min;
        }
        count = 0;
    }

    // Clustering time
    gettimeofday(&stop, 0);
    long sec = stop.tv_sec - start.tv_sec;
    long mic_sec = stop.tv_usec - start.tv_usec;
    double clustering_time = sec + mic_sec*1e-6;

    if(complete == 1)
    {
        for(int m=0; m<number_of_clusters; m++)
        {
            fprintf(output_file_ptr, "CLUSTER-%d ", m+1);
            fprintf(output_file_ptr, "{centroid: ");
            for(int i=0; i<dimension; i++)
                fprintf(output_file_ptr, "%f ", centroid_index[m][i]);
            for(int i=0; i<j[m]; i++)
                fprintf(output_file_ptr, ", %d", radius[m][i]);
            fprintf(output_file_ptr, "}\n\n");
        }
    }
    else
    {
        for(int m=0; m<number_of_clusters; m++)
        {
            fprintf(output_file_ptr, "CLUSTER-%d ", m+1);
            fprintf(output_file_ptr, "{size: %d, ", j[m]);
            fprintf(output_file_ptr, "centroid: ");
            for(int i=0; i<dimension; i++)
                fprintf(output_file_ptr, "%f ", centroid_index[m][i]);
            fprintf(output_file_ptr, "}\n\n");
        }
        fprintf(output_file_ptr, "clustering time: %f\n\n", clustering_time);
    }

    if(silhouette == 1)
    {
        // Silhouette
        long dist;
        int min_ind[number_of_clusters];
        for(int i=0; i<number_of_clusters; i++)     // which is the closest centroid to each centroid
        {
            min = INT_MAX;
            for(int j=0; j<number_of_clusters; j++)
            {
                if(i != j)
                {
                    dist = 0;
                    for(int d=1; d<=dimension; d++)
                    {
                        dist = dist + pow((centroid_index[i][d] - centroid_index[j][d]), 2);
                    }
                    dist = sqrt(dist);
                    if(dist < min)
                    {
                        min = dist;
                        min_ind[i] = j;
                    }
                }
            }
        }

        int clust, a_i, b;
        double sum = 0;
        int s[input_items_counter];
        for(int i=0; i<input_items_counter; i++)
        {
            clust = centroids[i];
            for(int k=0; k<j[clust]-1; k++)
            {
                dist = 0;
                for(int d=1; d<=dimension; d++)
                {
                    dist = dist + pow((curves[i][d-1] - curves[radius[clust][k]][d-1]), 2);   // distance for every item
                }
                dist = sqrt(dist);
                sum = sum + dist;
            }
            if(j[clust] == 0)
                a_i = 1;
            else
                a_i = sum / j[clust];   // average
            sum = 0;
            clust = min_ind[clust];     // closest centroid
            for(int k=0; k<j[clust]-1; k++)
            {
                dist = 0;
                for(int d=1; d<=dimension; d++)
                {
                    dist = dist + pow((curves[i][d-1] - curves[radius[clust][k]][d-1]), 2);
                }
                dist = sqrt(dist);
                sum = sum + dist;
            }
            if(j[clust] == 0)
                b = 1;
            else
                b = sum / j[clust]; // average
            sum = 0;

            if(a_i > b)
                s[i] = (b / a_i) - 1;
            else if(b > a_i)
                s[i] = 1 - (a_i / b);
            else
                s[i] = 0;
        }
        float average_per_cluster[number_of_clusters];  // average for each cluster
        for(int i=0; i<number_of_clusters; i++)
            average_per_cluster[i] = 0;
        for(int i=0; i<input_items_counter; i++)
        {
            for(int m=0; m<number_of_clusters; m++)
            {
                if(centroids[i] == m)
                {
                    average_per_cluster[m] = average_per_cluster[m] + s[i];
                    continue;
                }
            }
        }

        fprintf(output_file_ptr, "Silhouette: [");
        for(int m=0; m<number_of_clusters; m++)
        {
            fprintf(output_file_ptr, "%f, ", average_per_cluster[m] / j[m]);
        }

        float total_average_sil = 0;
        for(int i=0; i<input_items_counter; i++)
            total_average_sil = total_average_sil + s[i];
        total_average_sil = total_average_sil / input_items_counter;
        fprintf(output_file_ptr, "%f]\n\n\n", total_average_sil);
    }

    // free memory
    for(int i=0; i<input_items_counter; i++)
    {
        free(h_p_result[i]);
    }
    free(h_p_result);
    for(int i=0; i<number_of_clusters; i++)
    {
        free(previous_centroid_index[i]);
        // free(radius[i]);
        free(h_q_result[i]);
    }
    free(previous_centroid_index);
    // free(radius);
    free(h_q_result);
}


