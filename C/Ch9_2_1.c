#include <stdio.h>

int main(){
    int sum = 0, grade[4] = {81, 93, 77, 59};
    
    printf("grade_1 = %d\n", grade[0]);
    printf("grade_2 = %d\n", grade[1]);
    printf("grade_3 = %d\n", grade[2]);
    printf("grade_4 = %d\n", grade[3]);

    for(int i = 0; i < sizeof(grade) / sizeof(int); i++){
        sum += grade[i];
    }
    printf("sum = %d\n", sum);

    return 0;
}