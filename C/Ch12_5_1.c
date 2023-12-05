#include <stdio.h>

int main(){
    char a = 0x2c; /*00101100*/
    printf("a = %3d(%#x)\n", a, a);

    a = ~a;
    printf("~a = %3d(%#x)\n", a & 0xff, a & 0xff);

    a = ~a;
    printf("~(~a) = %3d(%#x)\n", a, a);
    
    return 0;
}