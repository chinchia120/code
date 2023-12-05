#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

class Solution 
{
public:
    int majorityElement(vector<int>& nums) 
    {
        if(nums.size() == 1)
        {
            return nums[0];
        }
        
        sort(nums.begin(), nums.end());
        int cnt = 0;
        vector<vector<int>> list = {{nums[0], 0}};

        for(int i = 0; i < nums.size()-1; i++)
        {
            list[cnt][1] += 1;

            if(nums[i] != nums[i+1])
            {
                cnt++;
                list.push_back(vector<int>());
                list[cnt].push_back(nums[i+1]);
                list[cnt].push_back(0);
            }          
        }

        if(nums[nums.size()-2] == nums[nums.size()-1])
        {
            list[cnt][1] += 1;
        }

        int max_number = 0, max_count = 0;
        for(int i = 0; i < list.size(); i++)
        {
            if(list[i][1] > max_count)
            {
                max_number = list[i][0];
                max_count = list[i][1];
            }
        }

        return max_number;
    }
};

int main(int argc, char **argv)
{
    vector<int> nums = {3, 2, 3};
    Solution S;
    int ans = S.majorityElement(nums);

    cout << ans;

    return 0;
}