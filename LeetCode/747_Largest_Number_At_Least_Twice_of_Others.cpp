#include <iostream> 
#include <vector>
#include <algorithm>
using namespace std;

class Solution 
{
public:
    int dominantIndex(vector<int>& nums) 
    {
        vector<int> nums_copy = nums;
        sort(nums.begin(), nums.end());
        if(nums[nums.size()-2] == 0 || nums[nums.size()-1]/nums[nums.size()-2] >= 2) 
        {
            return distance(nums_copy.begin(), find(nums_copy.begin(), nums_copy.end(), nums[nums.size()-1]));
        }
        else
        {
            return -1;
        }
    }
};

int main(int argc, char **argv)
{
    vector<int> nums = {3, 6, 1, 0};
    Solution S;

    cout << S.dominantIndex(nums);

    return 0;
}