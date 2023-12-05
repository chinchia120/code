#include <iostream>
#include <vector>
using namespace std;

class Solution 
{
public:
    int climbStairs(int n) 
    {
        vector<int> vec = {0, 1, 2};
        for(int i = 3; i <= n; i++)
        {
            vec.push_back(vec[i-1] + vec[i-2]);
        }

        return vec[n];
    }
};

int main(int argc, char **argv)
{
    int n = 3;
    Solution S;
    int ans = S.climbStairs(n);

    cout << ans;

    return 0;
}