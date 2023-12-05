#include <iostream>
#include <stdio.h>
using namespace std;

class Demo
{
public:
    int a;
    int b;
    string str1;
    string str2;
    int plus_opt();
    int minus_opt();
    int combine_opt();
};

int Demo::plus_opt()
{
    return a + b;
}

int Demo::minus_opt()
{
    return a - b;
}

int Demo::combine_opt()
{
    return (str1 + str2);
}

int main(int argc, char **argv)
{   
    Demo t;
    t.a = 11;
    t.b = 22;
    t.str1 = "test1,";
    t.str2 = " test2";
    
    printf("a = %d, b = %d\n", t.a, t.b);
    printf("plus_opt = %d\n", t.plus_opt());
    printf("minus_opt = %d\n", t.minus_opt());
    printf("combine_opt = %d\n", t.combine_opt());

    return 0;
}