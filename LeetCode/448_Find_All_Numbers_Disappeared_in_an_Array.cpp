#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

class Solution 
{
public:
    vector<int> findDisappearedNumbers(vector<int>& nums) 
    {   
        int len = nums.size();

        sort(nums.begin(), nums.end());
        //show_vector(nums);

        nums = remove_same_nums_in_vector(nums);
        //show_vector(nums);
  
        vector<int> dis_nums;
        int len_rm_vec = nums.size(), cnt = 0;
        for(int i = 1; i <= len; i++)
        {
            if(cnt < len_rm_vec)
            {
                if(nums[cnt] == i)
                {
                    cnt ++;
                }
                else
                {
                    dis_nums.push_back(i);
                }
            }
            else
            {
                dis_nums.push_back(i);
            }
        }

        return dis_nums;
    }

    void show_vector(vector<int> vec)
    {
        for(int i = 0; i < vec.size(); i++)
        {
            cout << vec[i] << " ";
        }
        cout << endl;
    }

    vector<int> remove_same_nums_in_vector(vector<int> vec)
    {   
        vector<int> after_remove_vector = {vec[0]};
        int cnt = 0;
        for(int i = 0; i < vec.size(); i++)
        {
            if(after_remove_vector[cnt] != vec[i])
            {
                after_remove_vector.push_back(vec[i]);
                cnt++;
            }
        }

        return after_remove_vector;
    }
};

int main(int argc, char **argv)
{
    vector<int> nums = {1, 1};
    Solution S;
    S.show_vector(S.findDisappearedNumbers(nums));

    return 0;
}