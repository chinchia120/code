#include <iostream>
#include <algorithm>
#include <vector>
using namespace std;

class Solution 
{
public:
    int singleNumber(vector<int>& nums) 
    {
        sort(nums.begin(), nums.end());
        nums.push_back(99999);

        int flag = 0, i = 0;
        for(i = 0; i < nums.size()-1; i+=2)
        {
            if(nums[i] != nums[i+1])
            {   
                flag = 1;
                break;
            }
        }
        
        return nums[i];
    }
};


int main(int argc, char ** argv)
{
    vector<int> nums = {4, 1, 2, 1, 2};
    Solution S;
    int ans = S.singleNumber(nums);

    cout << ans;

    return 0;
}