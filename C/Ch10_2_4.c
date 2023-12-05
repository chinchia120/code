#include <stdio.h>
#define length 5
void maxElement(int *, int *);

int main(){
    int index;
    int data[length] = {81, 93, 77, 59, 69};

    for(int i = 0; i < length; i++){
        printf("[%d : %d]", i, data[i]);
    }
    printf("\n");

    maxElement(data, &index);
    printf("max element = [%d : %d]\n", index, data[index]);

    return 0;
}

void maxElement(int *data, int *index){
    int maxValue = 0;

    for(int i = 0; i < length; i++){
        if(data[i] > maxValue){
            maxValue = data[i];
            *index = i;
        }
    }
}