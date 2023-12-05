#include <iostream>
using namespace std;

class Solution 
{
public:
    bool isPowerOfFour(int n) 
    {
        if(n == 0)
        {
            return false;
        }

        while(n != 1)
        {
            if(n%4 != 0)
            {
                return false;
            }

            n /= 4;
        }
        return true;
    }
};

int main(int argc, char **argv)
{
    int n = 27;
    Solution S;

    cout << S.isPowerOfFour(n) << endl;

    return 0;
}