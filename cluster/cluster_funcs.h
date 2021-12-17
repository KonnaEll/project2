#ifndef CLUSTER_FUNCS

int min_distance(double**, int, double**, int, int);
int max_distance(double*, int);
int min_distance_index(double**, int, double**, int, int);
void classic_assign(int, int, int, double**, double**, int, int, FILE*);
void lsh_assign(char**, int, int, int, double**, double**, int, int, FILE*, int);
void cube_assign(char** ,int, int, int, double**, double**, int, int, FILE*,int ,int ,int);
void frechet_classic_assign(int, int, int, double**, double**, int, int, FILE*);
void frechet_lsh_assign(char**, int, int, int, double**, double**, int, int, FILE*, int);

#endif