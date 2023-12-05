#include <stdio.h>

int main(){
    int step = 20, upper = 300;
    float f;

    for (int i = 0; i <= upper; i += step){
        f = (9.0 * i) / 5.0 + 32.0;
        
        printf("%d %.2f\n", i, f);
    }
    
    return 0;
}