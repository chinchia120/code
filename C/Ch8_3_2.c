#include <stdio.h>

void swap(int *, int *);

int main(){
    int a = 10, b = 15;
    printf("%d %d\n", a, b);
    swap(&a, &b);
    printf("%d %d\n", a, b);

    return 0;
}

void swap(int *x, int *y){
    int tmp = *x;
    *x = *y;
    *y = tmp;
}