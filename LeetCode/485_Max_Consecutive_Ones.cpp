#include <iostream>
#include <vector>
using namespace std;

class Solution 
{
public:
    int findMaxConsecutiveOnes(vector<int>& nums) 
    {
        int max = 0, sum = 0;

        for(int i = 0; i < nums.size(); i++)
        {
            if(nums[i] == 1)
            {
                sum++;
            }
            else
            {
                if(sum > max)
                {
                    max = sum;
                    sum = 0;
                }
            }
        }

        if(sum > max) max = sum;

        return max;
    }
};

int main(int argc, char **argv)
{
    vector<int> nums = {1, 1, 0, 1, 1, 1};
    Solution S;

    cout << S.findMaxConsecutiveOnes(nums);

    return 0;
}