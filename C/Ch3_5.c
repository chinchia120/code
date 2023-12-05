#include <stdio.h>
#define PI 3.1415926

int main()
{
    double area;
    double r = 10.0;

    const double e = 2.71828182845;
    area = PI * r * r;

    printf("Area = %f\n",area);
    printf("Const e = %f\n",e);

    return 0;
}