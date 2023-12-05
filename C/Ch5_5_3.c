#include <stdio.h>

int main()
{
    int a = 17, b = 5;
    float r;

    printf("a = %d, b = %d\n", a, b);
    
    r = a / b;
    printf("r = a / b = %f\n", r);

    r = (float)a / (float)b;
    printf("r = (float)a / (float)b = %f", r);

    return 0;

}