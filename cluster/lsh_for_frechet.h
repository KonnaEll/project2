#ifndef LSH_FOR_FRECHET

struct Node
{
    int x;
    double y;
};

float distance_computation(float**, int, double**, double**, int, int);
double** grid_to_frechet(double**, double, int, int, float*);
void lsh_for_frechet(double**, int, char**, double**, int, char**, int, double**, double**, FILE*, int, int);

#endif