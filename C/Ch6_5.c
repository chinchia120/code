#include <stdio.h>

int main(){
    int a = 3, b = 5, c = 2;

    if(a > b && a > c){
        printf("a is maximum value\n");
    }else{
        if(b > c){
            printf("b is maximum value\n");
        }else{
            printf("c is maximum value\n");
        }
    }

    return 0;
}