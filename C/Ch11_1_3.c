#include <stdio.h>
#include <string.h>

struct quiz{
    int math;
    int english;
};


struct student{
    int stdID;
    char name[20];
    struct quiz grade;
};

int main(){
    struct student std1;
    struct student std2 = {9402, "AAA", 65, 88};
    
    std1.stdID = 9401;
    strcpy(std1.name, "BBB");
    std1.grade.math = 90;
    std1.grade.english = 77;

    printf("student_ID     = %d\n", std1.stdID);
    printf("student_name   = %s\n", std1.name);
    printf("total score    = %d\n", std1.grade.math + std1.grade.english);
    printf("average score  = %.2f\n", (std1.grade.math + std1.grade.english) / 2.0);
    printf("-----------------------\n");

    printf("student_ID     = %d\n", std2.stdID);
    printf("student_name   = %s\n", std2.name);
    printf("total score    = %d\n", std2.grade.math + std2.grade.english);
    printf("average score  = %.2f\n", (std2.grade.math + std2.grade.english) / 2.0);
    printf("-----------------------\n");

    return 0;
}
