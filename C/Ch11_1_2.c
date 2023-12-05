#include <stdio.h>
#include <string.h>

struct student{
    int stdID;
    char name[20];
    int math;
    int english;
};

int main(){
    struct student std1;
    struct student std2 = {9402, "AAA", 65, 88};
    struct student std3;
    
    std1.stdID = 9401;
    strcpy(std1.name, "BBB");
    std1.math = 90;
    std1.english = 77;

    std3 = std2;

    printf("student_ID     = %d\n", std1.stdID);
    printf("student_name   = %s\n", std1.name);
    printf("total score    = %d\n", std1.math + std1.english);
    printf("average score  = %.2f\n", (std1.math + std1.english) / 2.0);
    printf("-----------------------\n");

    printf("student_ID     = %d\n", std2.stdID);
    printf("student_name   = %s\n", std2.name);
    printf("total score    = %d\n", std2.math + std2.english);
    printf("average score  = %.2f\n", (std2.math + std2.english) / 2.0);
    printf("-----------------------\n");

    printf("student_ID     = %d\n", std3.stdID);
    printf("student_name   = %s\n", std3.name);
    printf("total score    = %d\n", std3.math + std3.english);
    printf("average score  = %.2f\n", (std3.math + std3.english) / 2.0);
    printf("-----------------------\n");

    return 0;
}
