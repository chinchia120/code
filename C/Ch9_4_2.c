#include <stdio.h>
#define class 3
#define student 5

void maxScore(int [][student]);

int main(){
    int score[class][student] = {{74, 56, 33, 65, 89}, {37, 68, 44, 78, 92}, {33, 83, 77, 66, 88}};
    
    maxScore(score);

    return 0;
}

void maxScore(int score[][student]){
    int maxValue = 0, index_class, index_student;

    for(int i = 0; i < class; i++){
        for(int j = 0; j < student; j++){
            if(score[i][j] > maxValue){
                maxValue = score[i][j];
                index_class = i;
                index_student = j;
            }
        }
    }
    printf("class = %d, student = %d, grade = %d\n", index_class, index_student, maxValue);
}