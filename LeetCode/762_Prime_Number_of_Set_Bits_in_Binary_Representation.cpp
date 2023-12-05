#include <iostream>
using namespace std;

class Solution
{
public:
    int countPrimeSetBits(int left, int right)
    {   
        int cnt = 0;
        for(int i = left; i <= right; i++)
        {   
            int sum = 0, num = i;
            while(true)
            {   
                if(num % 2 == 1) sum++;
                num /= 2;
                if(num < 2)
                {
                    if(num == 1) sum++;
                    //cout << sum << endl;
                    break;
                }
            }
            if(sum == 2 || sum == 3 || sum == 5 || sum == 7 || sum == 11 || sum == 13 || sum == 17 || sum == 19) cnt++;
        }
        return cnt;
    }
};

int main(int argc, char **argv)
{
    int left = 6, right = 10;
    Solution S;

    cout << S.countPrimeSetBits(left, right);

    return 0;
}