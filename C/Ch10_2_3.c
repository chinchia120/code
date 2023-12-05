#include <stdio.h>

int main(){
    int var = 55, var1;
    int *ptr = NULL;

    ptr = &var;
    var1 = *ptr;

    printf("value of var = %d, address of var = %p\n", var, &var);
    printf("value of var1 = %d, address of var1 = %p\n", var1, &var1);
    printf("value of ptr = %p, address of ptr = %p\n", ptr, &ptr);
    printf("value of *ptr = %d\n", *ptr);

    return 0;
}