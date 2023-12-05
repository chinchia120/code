#include <stdio.h>
#define MAXSIZE 3

int main()
{
    struct test
    {
        int midterm;
        int final;
    };

    struct test student[MAXSIZE];

    for (int i = 0; i < MAXSIZE; i++)
    {
        printf("student_ID = %d\n", i + 1);
        printf("enter midterm test score ==> ");
        scanf("%d", &student[i].midterm);
        printf("enter  final  test score ==> ");
        scanf("%d", &student[i].final);
        printf("average score = %.2f\n", (student[i].midterm + student[i].final) / 2.0);
        printf("------------------------------------\n");
    }

    return 0;
}