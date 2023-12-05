#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

class Solution
{
public:
    vector<int> findErrorNums(vector<int>& nums)
    {
        sort(nums.begin(), nums.end());
        
        vector<int> ans;
        for(int i = 0; i < nums.size()-1; i++)
        {   
            if(nums[i] == nums[i+1]) ans.push_back(nums[i]);
        }
        
        vector<int> nums_rm = rm_same_vector(nums);
        
        int index = 0;
        for(int i = 1; i < nums.size()+1; i++)
        {   
            for(int j = index; j < nums_rm.size(); j++)
            {   
                //cout << i << " " << nums_rm[j] << endl;
                if(i < nums_rm[j])
                {
                    ans.push_back(i);
                    i++;
                }

                if(i == nums_rm[j])
                {   
                    index = j;
                    break; 
                }
            }
        }

        if(ans[ans.size()-1] != nums.size() && nums_rm[nums_rm.size()-1] != nums.size()) ans.push_back(nums.size());

        return ans;
    }

    void show_vector(vector<int> vec)
    {
        for(int i = 0; i < vec.size(); i++) cout << vec[i] << " ";
        cout << endl;
    }

    vector<int> rm_same_vector(vector<int> vec)
    {
        vector<int> rm_vec = {vec[0]};
        for(int i = 1; i < vec.size(); i++)
        {
            if(vec[i] != rm_vec[rm_vec.size()-1]) rm_vec.push_back(vec[i]);
        }
        return rm_vec;
    }
};

int main(int argc, char **argv)
{
    vector<int> nums = {3, 2, 2};
    Solution S;

    S.show_vector(S.findErrorNums(nums));

    return 0;
}
