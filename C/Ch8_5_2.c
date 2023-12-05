#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(){
    int rlt, num, times, playAgain = 1;
    srand(time(NULL));

    while(playAgain){
        rlt = rand()%100;
        times = 0;

        while(rlt != num){
            printf("enter a number between 0 to 99, ");
            scanf("%d", &num);
            times++;
            if(rlt == num){
                printf("get the answer %d in %d times\n", num, times);
                printf("do you want to play again ? ");
                scanf("%d", &playAgain);
            }else{
                if(rlt > num){
                    printf("bigger than %d\n", num);
                }else{
                    printf("smaller than %d\n", num);
                }
            }
        }
    }
    return 0;
}