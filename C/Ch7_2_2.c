#include <stdio.h>

int main(){
    int n, sum = 0;

    printf("enter a maximum number, ");
    scanf("%d", &n);

    for(int i = n; i > 0; i--){
        sum += i;
    }
    printf("sum = %d\n", sum);

    return 0;
}