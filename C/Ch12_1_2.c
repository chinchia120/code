#include <stdio.h>
#define FORMAT "number = %d\n"
#define MSG "end of program\n"
#define ONE 1
#define TWO ONE + ONE

int main()
{
    printf(FORMAT, ONE);
    printf(FORMAT, TWO);
    printf(MSG);

    return 0;
}