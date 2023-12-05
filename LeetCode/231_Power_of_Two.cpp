#include <iostream>
using namespace std;

class Solution 
{
public:
    bool isPowerOfTwo(int n) 
    {   
        if(n == 0)
        {
            return false;
        }

        if(n == 1)
        {
            return true;
        }

        while(n != 1)
        {
            int decimal = n % 2;

            if(decimal != 0)
            {   
                return false;
            }

            n /= 2;
        }

        return true;
    }
};

int main(int argc, char **argv)
{
    int n = 17;
    Solution S;
    cout << S.isPowerOfTwo(n);

    return 0;
}