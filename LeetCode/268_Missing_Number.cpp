#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

class Solution
{
public:
    int missingNumber(vector<int>& nums) 
    {
        sort(nums.begin(), nums.end());
        nums.push_back(nums[nums.size()-1]+2);

        if(nums[0] != 0) return 0;

        for(int i = 0; i < nums.size()-1; i++)
        {
            if(nums[i]+1 != nums[i+1]) return nums[i]+1;
        }

        return 0;
    }
};

int main(int argc, char **argv)
{
    vector<int> nums = {1, 2};
    Solution S;
    
    cout << S.missingNumber(nums) << endl;

    return 0;
}