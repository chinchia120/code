#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

class Solution
{
private:
    vector<int> remove_same(vector<int> vec)
    {   
        vector<int> rm_vec = {vec[0]};
        int cnt = 0;
        for(int i = 1; i < vec.size(); i++)
        {
            if(vec[i] != rm_vec[cnt])
            {
                rm_vec.push_back(vec[i]);
                cnt++;
            }
        }
        //show_vector(rm_vec);

        return rm_vec;
    }

    void show_vector(vector<int> vec)
    {
        for(int i = 0; i < vec.size(); i++)
        {
            cout << vec[i] << " ";
        }
        cout << endl;
    }

    vector<int> sort_vector(vector<int> vec)
    {
        for(int i = 0; i < vec.size(); i++)
        {
            for(int j = i+1; j < vec.size(); j++)
            {
                if(vec[i] > vec[j])
                {
                    int tmp = vec[j];
                    vec[j] = vec[i];
                    vec[i] = tmp;
                }
            }
        }
        //show_vector(vec);

        return vec;
    }

public:
    int longestConsecutive(vector<int>& nums) 
    {   
        if(nums.size() == 0)
        {
            return 0;
        }

        sort(nums.begin(), nums.end());
        //nums = sort_vector(nums);
        nums = remove_same(nums);

        if(nums.size() == 1)
        {
            return 1;
        }

        nums.push_back(nums[nums.size()-1]+10);
        int cnt = 0, tmp = 0;
        for(int i = 0; i < nums.size()-1; i++)
        {   
            cnt++;
            //cout << nums[i] << " " << nums[i+1] << endl;
            //cout << cnt << " " << tmp << endl;
            if(nums[i+1] != nums[i]+1)
            {
                if(cnt >= tmp)
                {
                    tmp = cnt;
                }
                cnt = 0;
            }
        }

        return tmp;
    }
};

int main(int argc, char **argv)
{
    vector<int> nums = {-7,-1,3,-9,-4,7,-3,2,4,9,4,-9,8,-7,5,-1,-7};
    Solution S;
    int ans = S.longestConsecutive(nums);

    cout << ans;

    return 0;
}