#include <stdio.h>
#define Length 6

int main(){
    int *ptr, data[Length] = {11, 93, 45, 27, -40, 80};

    for(int i = 0; i < Length; i++){
        ptr = &data[i];
        printf("data[%d] = %3d(%p)\n", i, *ptr, ptr);
    }

    ptr = data;
    printf("first element : %d(%p)\n", *ptr, &ptr);

    ptr = &data[Length - 1];
    printf("last element : %d(%p)\n", *ptr, &ptr);

    return 0;
}