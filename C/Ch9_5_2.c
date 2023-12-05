#include <stdio.h>
#include <string.h>

int main(){
    char str[] = "This is a pan.", str1[20], str2[20], str3[] = "Hi!";

    printf("str = %s\n", str);
    printf("length of str = %d\n",strlen(str));
    
    strcpy(str1, "This is a book.");
    strcpy(str2, str);

    printf("str1 = %s\n", str1);
    printf("str2 = %s\n", str2);

    strcat(str3, str1);
    printf("str3 = %s\n", str3);

    if(strcmp(str1, str2) > 0){
        printf("str1 is bigger.\n");
    }else{
        printf("str2 is bigger.\n");
    }

    return 0;
}