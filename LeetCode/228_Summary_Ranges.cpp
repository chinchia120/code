#include <iostream>
#include <vector>
#include <string>
using namespace std;

class Solution 
{
public:
    void show_vector(vector<string> vec)
    {
        for(int i = 0; i < vec.size(); i++)
        {
            cout << vec[i] << " ";
        }
        cout << endl;
    }

    vector<string> summaryRanges(vector<int>& nums) 
    {   
        vector<string> range;
        if(nums.size() == 1)
        {   
            range.push_back(to_string(nums[0]));

            return range;
        }
        else if(nums.size() == 0)
        {
            return range;
        }

        vector<int> check;
        nums.push_back(nums[nums.size()-1]);
        for(int i = 0; i < nums.size()-1; i++)
        {   
            check.push_back(nums[i]);
            if(nums[i] != nums[i+1]-1)
            {
                if(check.size() == 1)
                {   
                    range.push_back(to_string(check[0]));
                }   
                else
                {   
                    range.push_back(to_string(check[0]) + "->" + to_string(check[check.size()-1]));
                }
                check.clear();
            }
        }
        return range;
    }
};

int main(int argc, char **argv)
{
    vector<int> nums = {0, 1};
    Solution S;
    vector<string> ans = S.summaryRanges(nums);

    S.show_vector(ans);

    return 0;
}