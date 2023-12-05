#include <stdio.h>

int main(){
    int answer = 38, guess;

    while(1){
        printf("guess a number, ");
        scanf("%d", &guess);

        if(answer == guess){
            printf("got the answer %d\n", answer);
            break;
        }else{
            if(answer < guess){
                printf("smaller than %d\n", guess);
            }else{
                printf("bigger than %d\n", guess);
            }
        }
    }
    
    return 0;
}