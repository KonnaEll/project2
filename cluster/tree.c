#include <stdio.h>
#include <stdlib.h>
#include "tree.h"

struct Tree_Node* New_node(double* curve)
{
    struct Tree_Node* node = (struct Tree_Node*)malloc(sizeof(struct Tree_Node));
    node->curve = curve;
    node->left_child = NULL;
    node->right_child = NULL;
    return node;
}

struct Queue* New_Queue(int size) // create new queue
{
    struct Queue* queue = (struct Queue*)malloc(sizeof(struct Queue));

    queue->front = -1;
    queue->rear = -1;
    queue->size = size;

    queue->array = (struct Tree_Node**)malloc(sizeof(struct Tree_Node*) * queue->size);

    for(int i=0; i<size; i++)
        queue->array[i] = NULL;

    return queue;
}

int Empty(struct Queue* queue)
{
    return queue->front == -1;
}
 
int Full(struct Queue* queue)
{
    return queue->rear == queue->size - 1;
}
 
void Enqueue(struct Tree_Node* node, struct Queue* queue)
{
    if(Full(queue) == 1)
        return;
 
    queue->rear++;
    queue->array[queue->rear] = node;
 
    if(Empty(queue) == 1)
        queue->front++;
}

struct Tree_Node* Dequeue(struct Queue* queue)
{
    if(Empty(queue) == 1)
        return NULL;

    struct Tree_Node* temp = queue->array[queue->front];

    if (queue->front == queue->rear)
    {
        queue->front = -1;
        queue->rear = -1;
    }
    else
        queue->front++;

    return temp;
}

int Two_children(struct Tree_Node* temp)
{
    return temp && temp->left_child && temp->right_child;
}

void Insert_node(struct Tree_Node** root, double* curve, struct Queue* queue)
{
    struct Tree_Node* node = New_node(curve);

    if (*root == NULL)
        *root = node;

    else
    {
        struct Tree_Node* front = queue->array[queue->front];

        if(front->left_child == NULL)
            front->left_child = node;
        else if(front->right_child == NULL)
            front->right_child = node;

        if(Two_children(front))
            Dequeue(queue);
    }

    Enqueue(node, queue);
}

void levelOrder(struct Tree_Node* root)
{
    struct Queue* queue = New_Queue(100);

    Enqueue(root, queue);
    while (!Empty(queue))
    {
        struct Tree_Node* temp = Dequeue(queue);

        for(int i=0; i<729; i++)
            printf("%f\n", temp->curve[i]);

        if (temp->left_child)
            Enqueue(temp->left_child, queue);

        if (temp->right_child)
            Enqueue(temp->right_child, queue);
    }
}

double* Mean_curve(double* left_curve, double* right_curve, int dimension)
{
    double* mean = malloc(sizeof(double) * dimension);
    for(int i=0; i<dimension; i++)
        mean[i] = (left_curve[i] + right_curve[i]) / 2;
        
    return mean;
}


double* post_order_traversal(struct Tree_Node* node, int dimension)
{
    double* left_curve = malloc(sizeof(double) * dimension);
    double* right_curve = malloc(sizeof(double) * dimension);
    if(node->left_child == NULL && node->right_child == NULL)
        return node->curve;
    else
    {
        left_curve = post_order_traversal(node->left_child, dimension);
        if(node->right_child != NULL)
            right_curve = post_order_traversal(node->right_child, dimension);
        else
        {
            for(int i=0; i<dimension; i++)
                right_curve[i] = 0;
        }
    }

    return Mean_curve(left_curve, right_curve, dimension);
}