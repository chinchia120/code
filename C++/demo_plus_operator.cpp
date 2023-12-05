#include <iostream>
using namespace std;

int main(int argc, char **argv)
{   
    /* ----- initial value ----- */
    int i = 0, ans = 0;

    /* ----- plus operator ----- */
    i = 5;
    ans = i++;
    cout << "i = " << i << ", ans = "<< ans << endl; 

    i = 5;
    ans = ++i;
    cout << "i = " << i << ", ans = "<< ans << endl; 

    return 0;
}