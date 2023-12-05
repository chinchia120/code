#include <iostream>
using namespace std;

class Solution 
{
public:
    int arrangeCoins(int n) 
    {   
        int cnt = 0, plus = 1;
        while(n > 0)
        {   
            if(n > plus || n == plus) cnt++;
            n -= plus;
            plus ++;
        }
        return cnt;
    }
};

int main(int argc, char **argv)
{
    int n = 5;
    Solution S;
    
    cout << S.arrangeCoins(n);
    return 0;
}