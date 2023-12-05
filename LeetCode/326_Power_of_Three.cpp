#include <iostream>
using namespace std;

class Solution 
{
public:
    bool isPowerOfThree(int n) 
    {
        if(n == 0)
        {
            return false;
        }

        while(n != 1)
        {
            if(n%3 != 0)
            {
                return false;
            }

            n /= 3;
        }
        return true;
    }
};

int main(int argc, char **argv)
{
    int n = 27;
    Solution S;

    cout << S.isPowerOfThree(n) << endl;

    return 0;
}