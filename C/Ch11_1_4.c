#include <stdio.h>
#include <stdlib.h>

struct point{
    int x;
    int y;
};

struct point setXY(int, int);
struct point offset(struct point, int);

int main(){
    struct point p;
    struct point p1;
    
    p = setXY(150, 200);
    p1 = offset(p, 50);

    printf("original point = (%d, %d)\n", p.x, p.y);
    printf("  shift  point = (%d, %d)\n", p1.x, p1.y);

    return 0;
}

struct point setXY(int x, int y){
    struct point tmp;
    
    tmp.x= x;
    tmp.y = y;

    return tmp;
}

struct point offset(struct point p, int len){
    p.x += len;
    p.y += len;

    return p;
}