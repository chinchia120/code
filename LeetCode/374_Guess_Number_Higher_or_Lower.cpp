#include <iostream>
using namespace std;

class Solution 
{
public:
    int guessNumber(int n) 
    {   
        int min = 1, max = n;
        for(int i = min; i <= max; i++)
        {
            int repeat = guess(i);
            if(repeat == 0) return i;
            if(repeat == -1) max = i;
            if(repeat == 1) min = i;
        }
        return 0;
    }

    int guess(int num)
    {
        int pick = 6;
        if(num == pick) return 0;
        if(num > pick) return -1;
        if(num < pick) return 1;
    }
};

int main(int argc, char **argv)
{
    int n = 10;
    Solution S;

    cout << S.guessNumber(n);

    return 0;
}