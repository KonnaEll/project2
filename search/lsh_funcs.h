#ifndef LSH_FUNCS

struct Hash_Node
{
    char* name;
    int item;
    int ID;
    struct Hash_Node* next;
};

float h_function(double**, int, int);


#endif
