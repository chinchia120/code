#include <iostream>
#include <vector>
using namespace std;

class Solution 
{
public:
    int searchInsert(vector<int>& nums, int target) 
    {   
        if(nums[0] > target)
        {
            return 0;
        }

        if(nums.size() == 1)
        {
            if(nums[0] == target)
            {
                return 0;
            }
            else if(nums[0] > target)
            {
                return 0;
            }
            else
            {
                return 1;
            }
        }

        int i;
        for(i = 0; i < nums.size()-1; i++)
        {
            if(nums[i] == target)
            {
                return i;
            }
            
            if(nums[i] < target && nums[i+1] > target)
            {
                return i+1;
            }
        }

        if(nums[i] == target)
        {
            return i;
        }
        else
        {
            return i+1;
        }
    }
};

int main(int argc, char **argv)
{
    vector<int> demo = {1, 3, 5, 6};
    int target = 5;
    Solution S;
    int ans = S.searchInsert(demo, target);
    
    cout << ans;

    return 0;
}