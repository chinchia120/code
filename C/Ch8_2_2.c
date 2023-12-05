#include <stdio.h>

void sum(int, int);

int main(){

    sum(1, 5);
    sum(2, 7);

    return 0;
}

void sum(int start, int end){
    int sum = 0;

    for(int i = start; i <= end; i++){
        sum += i;
    }

    printf("add %d to %d equals to %d\n", start, end, sum);
}