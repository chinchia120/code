#include <stdio.h>

int main(){
    int *ptr = NULL;

    printf("value of ptr = %p, address of ptr = %p\n", ptr, &ptr);

    if(ptr == NULL){
        printf("ptr doesn't have initial value\n");
    }

    return 0;
}