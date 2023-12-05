#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

class Solution 
{
public:
    int thirdMax(vector<int>& nums) 
    {
        sort(nums.begin(), nums.end());
        nums = remove_same_nums_in_vector(nums);

        if(nums.size() < 3)
        {
            return nums[nums.size()-1];
        }

        return nums[nums.size()-3];
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
        //show_vector(after_remove_vector);

        return after_remove_vector;
    }

    void show_vector(vector<int> vec)
    {
        for(int i = 0; i < vec.size(); i++)
        {
            cout << vec[i] << " ";
        }
        cout << endl;
    }
};

int main(int argc, char **argv)
{
    vector<int> nums = {10, 9, 8, 8, 5};
    Solution S;

    cout << S.thirdMax(nums);

    return 0;
}