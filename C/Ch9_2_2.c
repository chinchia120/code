#include <stdio.h>

int main(){
    int sum = 0, score[] = {23, 32, 16, 22};

    for(int i = 0; i < sizeof(score) / sizeof(int); i++){
        sum += score[i];
    }
    printf("sum = %d\n", sum);
    printf("average = %d\n", sum/(sizeof(score) / sizeof(int)));
}