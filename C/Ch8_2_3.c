#include <stdio.h>

void printTriangle(int);
void isValid(double);
void C2F(double);

int main(){
    int n;
    double C;

    scanf("%d %lf", &n, &C);

    printTriangle(n);
    isValid(C);
    C2F(C);

    return 0;
}

void printTriangle(int row){
    for(int i = 1; i <= row; i++){
        for(int j = 1; j <= i; j++){
            printf("*");

            if(j == i){
                printf("\n");
            }
        }
    }
}

void isValid(double check){
    if(check >= 0 && check <= 200){
        printf("PASS\n");
    }else{
        printf("FAIL\n");
    }
}

void C2F(double C){
    printf("%.2lf\n", (9.0 * C) / 5.0 + 32.0);
}