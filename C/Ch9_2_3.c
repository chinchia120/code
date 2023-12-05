#include <stdio.h>
#define Length 4

int main(){
    double sum = 0, average, sale[Length] = {145.6, 178.9, 197.3, 156.7};

    for(int i = 0; i < Length; i++){
        sum += sale[i];
    }
    average = sum / Length;

    printf("sum = %lf\n", sum);
    printf("average = %lf\n", average);

    return 0;
}