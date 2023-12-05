#include <stdio.h>

int main(){
    int var = 100, var1 = 50;
    int *ptr = &var, *ptr1, *ptr2 = &var1;
    ptr1 = ptr;

    printf(" var  = %3d(%p)\n", var, &var);
    printf("*ptr  = %3d(%p)\n", *ptr, ptr);
    printf("*ptr1 = %3d(%p)\n", *ptr1, ptr1);
    printf("*ptr2 = %3d(%p)\n", *ptr2, ptr2);

    return 0;
}