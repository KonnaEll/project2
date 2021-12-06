#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lsh_funcs.h"
#include "cube_funcs.h"

static unsigned concatenate(unsigned x, unsigned y) {
    unsigned pow = 10;
    while(y >= pow)
        pow *= 10;                      //concat 2 unsigned numbers

    
    return x * pow + y;        
}

int * h_p(int** f, int k,int input_items_counter){
    int temp, a;
    int flag=0;
    int temp1=0;
    int *p = malloc(sizeof(int)*input_items_counter);
    for(long int i=0;i<input_items_counter;i++){
        for(int j=0;j<k;j++){
            if(flag==1){
                temp1++;
                a=concatenate(temp,f[i][j]);
                temp=a;
                if (temp1==k-1){
                    p[i]=temp;
                    temp1++;
                    flag=0;                     //concat d' (k) numbers 
                }                               // put them in table
            }else{
                temp=f[i][j];
            }
           
            if(temp1==k ){
                temp1=0;
            }else{
                flag=1;
            }
         }
     }

    return p;                                         //return table
}

int con(int x, int k){
    int digit;
    int i=0;
    int arr[k];
    int mynum=0;
    int temp=0;
    while (x!=0){
        digit=x % 10;
        arr[i]=digit;
        i++;                            
        x=x/10;
    }                                                   // converts a number from binary to decimal
    for(int j=0;j<i;j++)
    {
        temp=arr[j] * pow(2,j);
        mynum=mynum + temp;
    }

    return mynum;
}

int hammingDistance(int n1, int n2)
{
    int x = n1 ^ n2;
     int setBits = 0;

    while (x > 0)
    {                                   //calculate hamming distance between 2 integers
        setBits += x & 1;
        x >>= 1;
    }
 
    return setBits;
}