#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <limits.h>
#include "lsh_funcs.h"
#include "cube_funcs.h"
#include "cluster_funcs.h"

int main(int argc, char* argv[])
{
	// parameters check
    char* input_file = malloc(sizeof(char*) + 1);
    char* configuration_file = malloc(sizeof(char*) + 1);
    char* output_file = malloc(sizeof(char*) + 1);
    char* update = malloc(sizeof(char*) + 1);
    char* assignment = malloc(sizeof(char*) + 1);
    int complete = 0;
    int silhouette = 0;
	if(argc == 13 || argc == 12 || argc == 11)   // without complete
	{
        for(int i=1; i<=argc-1; i=i+2)
        {
            if(strcmp(argv[i], "-i") == 0)
                input_file = argv[i+1];
            else if(strcmp(argv[i], "-c") == 0)
                configuration_file = argv[i+1];
            else if(strcmp(argv[i], "-o") == 0)
                output_file = argv[i+1];
            else if(strcmp(argv[i], "-update") == 0)
                update = argv[i+1];
            else if(strcmp(argv[i], "-assignment") == 0)
                assignment = argv[i+1];
            else if(strcmp(argv[i], "-complete") == 0)
                complete = 1;
            else if(strcmp(argv[i], "-silhouette") == 0)
                silhouette = 1;                
        }
    }
    else
    {
        printf("Wrong input. Try again!\n");
    }

    printf("%s %s %s %s %s\n", input_file, configuration_file, output_file, update, assignment);


    FILE *input_file_ptr;
    input_file_ptr = fopen(input_file, "r");    // open input file
    if(input_file_ptr == NULL)
    {
        perror("Error\n");
        exit(1);
    }

    char c;
    int input_items_counter = 0;
    for (c = getc(input_file_ptr); c != EOF; c = getc(input_file_ptr))  // count the items of the file input(axis x)
    {
        if (c == '\n')
            input_items_counter = input_items_counter + 1;
    }
    rewind(input_file_ptr);

    int dimension = 0;
    while((c = fgetc(input_file_ptr)) != '\n')  // calculate the dimension of all the vectors. One vector is enough.
    {
        if(c == '\t')
            dimension++;
    }
    rewind(input_file_ptr);

    // Work for input file
    double** curves = malloc(sizeof(double*) * input_items_counter);    // array of the items of dataset
    for(int i=0; i<input_items_counter; i++)
        curves[i] = malloc(sizeof(double) * dimension);

    char** names = malloc(sizeof(char*) * input_items_counter);
    for(int i=0; i<input_items_counter; i++)
        names[i] = malloc(sizeof(char) * 20);

    int i = 0, j = 0, n = 0;
    char* x = malloc(sizeof(char*) + 1);
    fscanf(input_file_ptr, "%s", x);
    strcpy(names[n], x);
    while(fscanf(input_file_ptr, "%s", x) != EOF)     // fill the array with the dataset
    {
        curves[i][j] = atof(x);
        j++;
        if(j == dimension)
        {
            i++;
            j = 0;
            n++;
            if(fscanf(input_file_ptr, "%s", x) != EOF)
                strcpy(names[n], x);
        }
    }

    FILE *output_file_ptr;
    output_file_ptr = fopen(output_file, "w");    // open output file
    if(output_file_ptr == NULL)
    {
        perror("Error\n");
        exit(1);
    }

    FILE *conf_file_ptr;
    conf_file_ptr = fopen(configuration_file, "r");    // open configuration file
    if(conf_file_ptr == NULL)
    {
        perror("Error\n");
        exit(1);
    }

    int number_of_clusters;
    int number_of_vector_hash_tables;
    int number_of_vector_hash_functions;
    int max_number_M_hypercube;
    int number_of_hypercube_dimensions;
    int number_of_probes;
    i = 0;
    char* numb = malloc(sizeof(char*) + 1);
    while(fscanf(conf_file_ptr, "%s", numb) != EOF)     // save the contents of the file in variables
    {
        if(i == 3)
            number_of_clusters = atoi(numb);
        else if(i == 6)
            number_of_vector_hash_tables = atoi(numb);
        else if(i == 9)
            number_of_vector_hash_functions = atoi(numb);
        else if(i == 12)
            max_number_M_hypercube = atoi(numb);
        else if(i == 15)
            number_of_hypercube_dimensions = atoi(numb);
        else if(i == 18)
            number_of_probes = atoi(numb);
        i++;
    }
    printf("%d %d %d %d %d %d %d %d\n", complete, silhouette, number_of_clusters, number_of_vector_hash_tables, number_of_vector_hash_functions, number_of_probes, number_of_hypercube_dimensions, max_number_M_hypercube);

    // Initialize the centroids with K-means++ method
    srand(time(0));
    // array with number of clusters lines and dimension columns
    double** centroid_index = malloc(sizeof(double*) * number_of_clusters);
    for(int i=0; i<number_of_clusters; i++)
    {
        centroid_index[i] = malloc(sizeof(double) * (dimension + 1));
    }

    int item_id = rand() % input_items_counter;
    centroid_index[0][0] = item_id;
    for(int d=1; d<=dimension; d++)
    {
        centroid_index[0][d] = curves[item_id][d-1];   // put the dimensions
    }

    double* D = malloc(sizeof(double) * input_items_counter);
    int index;
	for(int i=1; i<number_of_clusters; i++)
    {
		for(int j=0; j<input_items_counter; j++)
        {
			D[j] = min_distance(curves, j, centroid_index, i, dimension);    // min distance from the centroids
		}

        index = max_distance(D, input_items_counter);    // max distance of the min distances above
        centroid_index[i][0] = index;

        for(int d=1; d<=dimension; d++)
        {
            centroid_index[i][d] = curves[index][d-1];
        }
	}
    // for(int i=0; i<number_of_clusters; i++)
    // {
    //     printf("\n\n");
    //     for(int d=0; d<=dimension; d++)
    //         printf("%f ", centroid_index[i][d]);
    // }

    if(strcmp(assignment, "Classic") == 0)
    {
        classic_assign(input_items_counter, dimension, number_of_clusters, curves, centroid_index, complete, silhouette, output_file_ptr);
    }
    else if(strcmp(assignment, "LSH") == 0)
    {
        lsh_assign(names, input_items_counter, dimension, number_of_clusters, curves, centroid_index, complete, silhouette, output_file_ptr, number_of_vector_hash_functions);
    }
    else if(strcmp(assignment, "Hypercube") == 0)
    {
        cube_assign(names, input_items_counter, dimension, number_of_clusters, curves, centroid_index, complete, silhouette, output_file_ptr, number_of_vector_hash_functions,number_of_probes,max_number_M_hypercube);
    }
    else if(strcmp(assignment, "LSH_Frechet") == 0)
    {

    }













    return 0;
}