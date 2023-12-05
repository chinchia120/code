#include <stdio.h>

void change(int *, int);

int main(){
    int data[] = {81, 13, 27, 39, 69};
    int length = sizeof(data)/sizeof(int);

    printf("before calling finction, ");
    for(int i = 0; i < length; i++){
        printf("[%d, %d] ", i, data[i]);
    }
    printf("\n");

    change(data, length);

    printf("after calling finction, ");
    for(int i = 0; i < length; i++){
        printf("[%d, %d] ", i, data[i]);
    }
    printf("\n");
    
    return 0;
}

void change(int data[], int length){
    int minValue = 100, index;

    for(int i = 0; i < length; i++){
        if(data[i] < minValue){
            minValue = data[i];
            index = i;
        }
    }
    
    data[index] = data[0];
    data[0] = minValue;   
}