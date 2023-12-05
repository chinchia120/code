#include <iostream>
using namespace std;

class Solution 
{
public:
    int mySqrt(int x)
    {   
        for(long long int i = 0; i <= x; i++)
        {   
            if(i*i == x)
            {   
                return i;
            }
            
            if(i*i < x && (i+1)*(i+1) > x)
            {   
                return i;
            }
        }

        return 0;
    }
};

int main(int argc, char **argv)
{
    int x = 0;
    Solution S;
    int ans = S.mySqrt(x);

    cout << ans;

    return 0;
}