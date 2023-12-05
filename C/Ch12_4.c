#include <stdio.h>

int main()
{
    char a = 0x3c; /*00111100*/
    char b = 0x0f; /*00001111*/
    char c = 0x03; /*00000011*/
    char d = 0x78; /*01111000*/
    char r;

    printf("a = %3d(%#.2x)\n", a, a);
    printf("b = %3d(%#.2x)\n", b, b);
    printf("c = %3d(%#.2x)\n", c, c);
    printf("d = %3d(%#.2x)\n", d, d);

    r = a & b;
    printf("a & b = %3d(%#.2x)\n", r, r);

    r = a | c;
    printf("a | c = %3d(%#.2x)\n", r, r);

    r = a ^ d;
    printf("a ^ d = %3d(%#.2x)\n", r, r);

    return 0;
}