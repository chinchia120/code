#include <stdio.h>
#include "Ch12_1_1.h"

int main()
{
    int r;

    do
    {
        printf("enter a radius in circle ==> ");
        scanf("%d", &r);

        if (r >= 0)
        {
            printf("circular area = %.3f\n", area(r));
        }
    } while (r >= 0);

    return 0;
}

double area(int r)
{
    return PI * r * r;
}