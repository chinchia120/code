#include <stdio.h>

void funcA();
void funcB();

int a, b = 2;

int main(){
    printf("initial value, a(G) = %d, b(G) = %d\n", a, b);
    funcA();
    printf("after calling funcA(), a(G) = %d, b(G) = %d\n", a, b);
    funcB();
    printf("after calling funcA(), a(G) = %d, b(G) = %d\n", a, b);
    
    return 0;
}

void funcA(){
    int a = 3;
    printf("in funcA, a(L) = %d, b(G) = %d\n", a, b);
    printf("a + b = %d\n",a + b);
}

void funcB(){
    a = 3;
    b = 4;
    printf("in funcB, a(G) = %d, b(G) = %d\n", a, b);
    printf("a + b = %d\n",a + b);
}