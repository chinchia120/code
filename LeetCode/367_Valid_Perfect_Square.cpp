#include <iostream>
using namespace std;

class Solution 
{
public:
    bool isPerfectSquare(int num) 
    {
        for(long long int i = 1; i <= num; i++)
        {
            if(i*i == num)
            {
                return true;
            }

            if(i*i > num)
            {
                return false;
            }
        }

        return false;
    }
};

int main(int argc, char **argv)
{
    int num = 17;
    Solution S;

    cout << S.isPerfectSquare(num);

    return 0;
}