#ifndef TREE

struct Tree_Node
{
    double* curve;
    struct Tree_Node* left_child;
    struct Tree_Node* right_child;
};

struct Queue
{
    int front;
    int rear;
    int size;
    struct Tree_Node** array;
};

struct Tree_Node* New_Node(double*);
// struct Tree_Node* Insert_right_Tree_Node(struct Tree_Node*, int);
// struct Tree_Node* Insert_left_Tree_Node(struct Tree_Node*, int);
struct Queue* New_Queue(int);
void Enqueue(struct Tree_Node*, struct Queue*);
struct Tree_Node* Dequeue(struct Queue*);
int Two_children(struct Tree_Node*);
void Insert_node(struct Tree_Node**, double*, struct Queue*);
void levelOrder(struct Tree_Node*);
int Empty(struct Queue* queue);
int Full(struct Queue* queue);

double* Mean_curve(double*, double*, int);
double* post_order_traversal(struct Tree_Node*, int);

#endif