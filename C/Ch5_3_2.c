#include <stdio.h>

int main()
{
    int c;
    double f;

    printf("enter degree of Celsius, ");
    scanf("%d", &c);

    f = (9.0 * c) / 5.0 + 32.0;
    printf("equal to %.2lf degree Fahrenheit.\n", f);

    return 0;
}