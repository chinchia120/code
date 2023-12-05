#include <iostream>
using namespace std;

int main(int argc, char **argv)
{   
    /* ----- initial value ----- */
    int i = 0;

    /* ----- for loop ----- */
    for(i = 0; i < 10; i++)
    {
        cout << i << " ";
    }
    cout << endl << "i = " << i << endl;

    for(i = 0; i < 10; ++i)
    {
        cout << i << " ";
    }
    cout << endl << "i = " << i << endl;
    
    return 0;
}