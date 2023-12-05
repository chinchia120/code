#include <stdio.h>

int main(){
    int id, sum = 0, grades[3][2] = {{74, 56}, {37, 68}, {33, 83}};
    double average;

    printf("enter the id, ");
    scanf("%d", &id);

    if(id >= 0 && id <= 2){
        for(int i = 0; i < 2; i++){
            sum += grades[id][i];
            printf("grade_%d = %d\n", i+1, grades[id][i]);
        }
        printf("id = %d, sum = %d\n", id, sum);

        average = sum / 2.0;
        printf("id = %d, average = %.2lf\n", id, average);
    }

    return 0;
}