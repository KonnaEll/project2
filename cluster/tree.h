#ifndef TREE

struct Node
{
    double* curve;
    struct Node* left_child;
    struct Node* right_child;
};

struct Queue
{
    int front;
    int rear;
    int size;
    struct Node** array;
};

struct Node* New_node(double*);
// struct Node* Insert_right_node(struct Node*, int);
// struct Node* Insert_left_node(struct Node*, int);
struct Queue* New_Queue(int);
void Enqueue(struct Node*, struct Queue*);
struct Node* Dequeue(struct Queue*);
int Two_children(struct Node*);
void Insert_node(struct Node**, double*, struct Queue*);
void levelOrder(struct Node*);
int Empty(struct Queue* queue);
int Full(struct Queue* queue);

double* Mean_curve(double*, double*, int);
double* post_order_traversal(struct Node*, int);

#endif