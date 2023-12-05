#include <stdio.h>

int main()
{
    int var1, var2;
    
    printf("enter operand_1, ");
    scanf("%d", &var1);

    printf("enter operand_2, ");
    scanf("%d", &var2);
   
    printf("after multiply, %d\n", var1 * var2);
    printf("after divided, %d\n", var1 / var2);

    return 0;
}