#include <stdio.h>

void printMsg(){
    printf("Welecome\n");
}

void sum2Ten(){
    int sum = 0;

    for(int i = 1; i < 11; i++){
        sum += i;
    }

    printf("sum = %d\n", sum);
}

int main(){
    printMsg();
    sum2Ten();

    return 0;
}