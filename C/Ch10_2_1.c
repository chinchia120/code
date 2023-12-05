#include <stdio.h>

int main(){
    int var = 100;
    int *ptr = &var;

    printf("valve of var = %d, address of var = %p\n", var, &var);
    printf("valve of ptr = %p, address of ptr = %p\n", ptr, &ptr);

    return 0;
}