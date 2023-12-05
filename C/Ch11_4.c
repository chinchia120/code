#include <stdio.h>

int main(){
    union number{
        char c;
        short value;
    };

    union number no;
    short ch, num;
    printf("size of the union = %d byte\n", sizeof(union number));

    do{
        printf("enter a number in Hexadecimal ==> ");
        scanf("%x", &num);

        if(num == -1){
            break;
        }

        no.value = num;

        printf("no.value = %#x(%d)\n", no.value, no.value);
        printf("last 8 bit = %d\n", no.value & 0x00ff);
        printf("no.c = %c(%d)\n", no.c, no.c);

        while((ch = getchar()) != '\n');

        printf("enter a char ==> ", ch);
        scanf("%c", &ch);
        no.c = (char)ch;

        printf("no.value = %#x(%d)\n", no.value, no.value);
        printf("last 8 bit = %d\n", no.value & 0x00ff);
        printf("no.c = %c(%d)\n", no.c, no.c);


    }while(1);

    return 0;
}