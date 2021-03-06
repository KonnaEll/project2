#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <limits.h>
#include "lsh_funcs.h"
#include "cube_funcs.h"
#include "lsh_for_vectors.h"
#include "cube_for_vectors.h"
#include "lsh_for_frechet.h"
#include "lsh_for_continuous.h"

int main(int argc, char* argv[])
{
    // parameters check
    char* input_file = malloc(sizeof(char*) + 1);
    char* query_file = malloc(sizeof(char*) + 1);
    char* output_file = malloc(sizeof(char*) + 1);
    char* algorithm = malloc(sizeof(char*) + 1);
    char* metric = malloc(sizeof(char*) + 1);
    int k = -1;
    int L = -1;
    int M = -1;
    int probes = -1;
    double delta;
    
    int inp_count = 0, que_count = 0, out_count = 0, alg_count = 0, met_count = 0;
    for(int i=1; i<=argc-2; i=i+2)
    {
        if(strcmp(argv[i], "-i") == 0)
        {
            input_file = argv[i+1];
            inp_count++;
        }
        else if(strcmp(argv[i], "-q") == 0)
        {
            query_file = argv[i+1];
            que_count++;
        }
        else if(strcmp(argv[i], "-o") == 0)
        {
            output_file = argv[i+1];
            out_count++;
        }
        else if(strcmp(argv[i], "-algorithm") == 0)
        {
            algorithm = argv[i+1];
            alg_count++;
        }
        else if(strcmp(argv[i], "-metric") == 0)
        {
            metric = argv[i+1];
            met_count++;
        }
        else if(strcmp(argv[i], "-k") == 0)
            k = atoi(argv[i+1]);
        else if(strcmp(argv[i], "-L") == 0)
            L = atoi(argv[i+1]);
        else if(strcmp(argv[i], "-M") == 0)
            M = atoi(argv[i+1]);
        else if(strcmp(argv[i], "-probes") == 0)
            probes = atoi(argv[i+1]);
        else if(strcmp(argv[i], "-delta") == 0)
            delta = atof(argv[i+1]);                                    
    }

    if(k == -1)
        k = 4;
    if(L == -1)
        L = 5;
    if(M == -1)
        M = 10;
    if(probes == -1)
        probes = 2;

    printf("\n");
    if(inp_count == 0)
    {
        printf("Give dataset path\n");
        scanf("%s", input_file);
    }
    if(que_count == 0)
    {
        printf("Give query file\n");
        scanf("%s", query_file);
    }
    if(out_count == 0)
    {
        printf("Give output file\n");
        scanf("%s", output_file);
    }
    if(alg_count == 0)
    {
        printf("Give algorithm\n");
        scanf("%s", algorithm);
    }
    if(met_count == 0)
    {
        printf("Give metric\n");
        scanf("%s", metric);
    }
    if((strcmp(algorithm, "Frechet") != 0) && (strcmp(metric, "continuous") == 0))
    {
        printf("Wrong input. No metric 'continuous' for LSH or Hypercube\n");
        exit(1);
    }

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
    while((c = fgetc(input_file_ptr)) != '\n')  // calculate the dimension of all the curves. One curve is enough.
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

    // for query file
    FILE *query_file_ptr;
    query_file_ptr = fopen(query_file, "r");    // open query file
    if(query_file_ptr == NULL)
    {
        perror("Error\n");
        exit(1);
    }
    int query_items_counter = 0;
    for (c = getc(query_file_ptr); c != EOF; c = getc(query_file_ptr))  // count the items of the file query(axis x)
    {
        if (c == '\n')
            query_items_counter = query_items_counter + 1;
    }
    rewind(query_file_ptr);

    // Work for query file
    double** query_curves = malloc(sizeof(double*) * query_items_counter);    // array of the items of dataset
    for(int i=0; i<query_items_counter; i++)
        query_curves[i] = malloc(sizeof(double) * (dimension + 1));

    char** query_names = malloc(sizeof(char*) * query_items_counter);
    for(int i=0; i<query_items_counter; i++)
        query_names[i] = malloc(sizeof(char) * 20);

    i = 0, j = 0, n = 0;
    char* y = malloc(sizeof(char*) + 1);
    fscanf(query_file_ptr, "%s", y);
    strcpy(query_names[n], y);
    while(fscanf(query_file_ptr, "%s", y) != EOF)     // fill the array with the dataset
    {
        query_curves[i][j] = atof(y);
        j++;
        if(j == dimension)
        {
            i++;
            j = 0;
            n++;
            if(fscanf(query_file_ptr, "%s", y) != EOF)
                strcpy(query_names[n], y);
        }
    }

    FILE *output_file_ptr;
    output_file_ptr = fopen(output_file, "w");    // open output file
    if(output_file_ptr == NULL)
    {
        perror("Error\n");
        exit(1);
    }

    srand(time(0));
    if(strcmp(algorithm, "LSH") == 0)
    {
        for(int n=0; n<L; n++)
        {
            lsh_for_vectors(input_items_counter, names, query_items_counter, query_names, dimension, curves, query_curves, output_file_ptr, k, n);
        }
    }
    else if(strcmp(algorithm, "Hypercube") == 0)
    {
        cube_for_vectors(input_items_counter, names, query_items_counter, query_names, dimension, curves, query_curves, output_file_ptr, k,probes,M);
    }
    else if(strcmp(algorithm, "Frechet") == 0)
    {
        if(strcmp(metric, "discrete") == 0)
        {
            double** vectors = malloc(sizeof(double*) * input_items_counter);    // array of the items of dataset
            for(int i=0; i<input_items_counter; i++)
                vectors[i] = malloc(sizeof(double) * (2*dimension + 1));

            double** query_vectors = malloc(sizeof(double*) * query_items_counter);    // array of the items of dataset
            for(int i=0; i<query_items_counter; i++)
                query_vectors[i] = malloc(sizeof(double) * (2*dimension + 1));

            float* t = malloc(sizeof(float) * dimension);
            for(int n=0; n<L; n++)
            {
                for(int i=0; i<dimension; i++)
                    t[i] = ((float)rand() / (float)(RAND_MAX)) * delta;

                vectors = grid_to_frechet(curves, delta, input_items_counter, dimension, t);
                query_vectors = grid_to_frechet(query_curves, delta, query_items_counter, dimension, t);

                lsh_for_frechet(vectors, input_items_counter, names, query_vectors, query_items_counter, query_names, dimension, curves, query_curves, output_file_ptr, k, n);
            }
        }
        else if(strcmp(metric, "continuous") == 0)
        {
            double** vectors = malloc(sizeof(double*) * input_items_counter);    // array of the items of dataset
            for(int i=0; i<input_items_counter; i++)
                vectors[i] = malloc(sizeof(double) * (dimension + 1));

            double** query_vectors = malloc(sizeof(double*) * query_items_counter);    // array of the items of dataset
            for(int i=0; i<query_items_counter; i++)
                query_vectors[i] = malloc(sizeof(double) * (dimension + 1));

            double** filtered = malloc(sizeof(double*) * input_items_counter);
            for(int i=0; i<input_items_counter; i++)
                filtered[i] = malloc(sizeof(double) * dimension);
            
            double** query_filtered = malloc(sizeof(double*) * query_items_counter);
            for(int i=0; i<query_items_counter; i++)
                query_filtered[i] = malloc(sizeof(double) * dimension);

            double** timeseries = malloc(sizeof(double*) * input_items_counter);    // array of the items of dataset
            for(int i=0; i<input_items_counter; i++)
                timeseries[i] = malloc(sizeof(double) * (dimension + 1));    

            double** query_timeseries = malloc(sizeof(double*) * query_items_counter);    // array of the items of dataset
            for(int i=0; i<query_items_counter; i++)
                query_timeseries[i] = malloc(sizeof(double) * (dimension + 1)); 


            filtered = filtering(curves, input_items_counter, dimension);
            query_filtered = filtering(query_curves, query_items_counter, dimension);

            vectors = grid_to_vector(filtered, delta, input_items_counter, dimension);
            query_vectors = grid_to_vector(query_filtered, delta, query_items_counter, dimension);
            
            timeseries = min_max_padding(vectors, input_items_counter, dimension);
            query_timeseries = min_max_padding(query_vectors, query_items_counter, dimension);

            // for(int j=0; j<dimension; j++)
            //     printf("%f ", vectors[50][j]);
            // printf("\n\n\n");
            // for(int j=0; j<dimension; j++)
            //         printf("%f ", timeseries[50][j]);
            // printf("\n");

            lsh_for_continuous(timeseries, input_items_counter, names, query_timeseries, query_items_counter, query_names, dimension, curves, query_curves, output_file_ptr, k);
        }
    }

    return 0;
}