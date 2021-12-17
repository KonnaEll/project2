#ifndef LSH_FOR_CONTINUOUS

double** filtering(double**, int, int);
double** grid_to_vector(double**, double, int, int);
double** min_max_padding(double**, int, int);
void lsh_for_continuous(double**, int, char**, double**, int, char**, int, double**, double**, FILE*, int);

#endif