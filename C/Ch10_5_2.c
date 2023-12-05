#include <stdio.h>

int main(){
    char str[15] = "This is a pen.", str1[15] = "";
    char *ptr, *ptr1;
    int i = 0;

    ptr = str;
    ptr1 = str1;
    

    while(*ptr != '\0'){
        *(ptr1 + i) = *ptr++;
        i++;
    }
    *(ptr + i) = '\0';

    printf("str  = %s\n", str);
    printf("str1 = %s\n", str1);
    printf("ptr1 = %s\n", ptr1);

    return 0;
}