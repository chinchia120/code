#include <stdio.h>

int main(){
    int data[] = {1, 2, 3, 4, 5}, i, j;
    int *ptr1 = &data[0], *ptr2 = &data[4];

    printf("data[0] = %d(%p)\n", data[0], &data[0]);
    printf("data[1] = %d(%p)\n", data[1], &data[1]);
    printf("data[2] = %d(%p)\n", data[2], &data[2]);
    printf("data[3] = %d(%p)\n", data[3], &data[3]);
    printf("data[4] = %d(%p)\n", data[4], &data[4]);
    printf("\n");

    printf("pointer operator : \n");

    printf("*ptr1 = %d(%p)\n", *ptr1, ptr1);
    printf("*ptr2 = %d(%p)\n", *ptr2, ptr2);
    printf("\n");

    ptr1++;
    printf("*ptr1++ = %d(%p)\n", *ptr1, ptr1);

    ptr1--;
    printf("*ptr1-- = %d(%p)\n", *ptr1, ptr1);

    ptr1 += 3;
    printf("*ptr1-- = %d(%p)\n", *ptr1, ptr1);

    ptr1 -= 2;
    printf("*ptr1-- = %d(%p)\n", *ptr1, ptr1);

    i = (int)(ptr2 - ptr1);
    j = (int)(ptr1 - ptr2);
    printf("ptr2 - ptr1 = %2d\n", i);
    printf("ptr1 - ptr2 = %2d\n", j);

    printf("\n");
    ptr1 = &data[0];
    for(int i = 0; i < 5; i++){
        printf("data[%d] = %d\t", i, *ptr1++);
    }
    printf("\n");

    return 0;
}