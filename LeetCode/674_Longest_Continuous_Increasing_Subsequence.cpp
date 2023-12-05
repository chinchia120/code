#include <iostream> 
#include <vector>
using namespace std;

class Solution 
{
public:
    int findLengthOfLCIS(vector<int>& nums) 
    {
        int len = 1, cnt = 1, tmp = nums[0];
        for(int i = 1; i < nums.size(); i++)
        {
            if(nums[i] > tmp) cnt++;
            else cnt = 1;
            tmp = nums[i];
            len = max(len, cnt);
        }

        return len;
    }
};

int main(int argc, char **argv)
{
    vector<int> nums = {1,3,5,4,7};
    Solution S;

    cout << S.findLengthOfLCIS(nums);

    return 0;
}