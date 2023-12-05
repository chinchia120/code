#include <stdio.h>

int main(){
    char str1[15] = "This is a pen.";
    char str2[15] = "hello! world";
    char *ptr1, *ptr2;

    printf("str1 = %s\n", str1);
    printf("str2 = %s\n", str2);
    ptr1 = str1;
    ptr2 = str2;
    printf("ptr1 = %s\n", ptr1);
    printf("ptr2 = %s\n", ptr2);

    return 0;
}