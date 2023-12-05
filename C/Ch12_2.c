#include <stdio.h>
#define TRUE 1
#define FALSE 0
#define SWAP(x, y){ \
    int _z; \
    _z = y; \
    y = x; \
    x = _z; \
}

#define AREA(r) (r*r*3.1415926)
#define ISEVEN(x) (x % 2 == 0) ? TRUE : FALSE

/*
#define ISEVEN(x){ \
    if(x % 2 == 0){ \
        return TRUE; \
    }else{ \
        return FALSE; \
    } \
}
*/

int main(){
    int a = 10, b = 5, r = 10;

    if(ISEVEN(a + b)){
        printf("a + b is EVEN\n");
    }else{
        printf("a + b is ODD\n");
    }

    printf("before swap a = %2d, b = %2d\n", a, b);

    SWAP(a, b);
    printf("after  swap a = %2d, b = %2d\n", a ,b);
    printf("radius = %d, circular area = %.3f\n", r, AREA(r));
    printf("radius = 25, circular area = %.3f\n", AREA(25));

    return 0;
}