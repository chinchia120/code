#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

class Solution 
{
public:
    int maximumProduct(vector<int>& nums) 
    {   
        if(nums.size() == 3) return nums[0] * nums[1] * nums[2];

        //show_vector(nums);
        sort(nums.begin(), nums.end());
        //show_vector(nums);
        int len = nums.size(), max_sum = nums[0] * nums[1] * nums[2];
        if(nums[0] * nums[1] * nums[len-1] > max_sum) max_sum = nums[0] * nums[1] * nums[len-1];
        if(nums[0] * nums[len-2] * nums[len-1] > max_sum) max_sum = nums[0] * nums[len-2] * nums[len-1];
        if(nums[len-3] * nums[len-2] * nums[len-1] > max_sum) max_sum = nums[len-3] * nums[len-2] * nums[len-1];

        return max_sum;
    }

    void show_vector(vector<int> vec)
    {
        for(int i = 0; i < vec.size(); i++) cout << vec[i] << " ";
        cout << endl;
    }
};

int main(int argc, char ** argv)
{
    //vector<int> nums = {-1, -2, -3};
    vector<int> nums = {1, 2, 3, 4};
    Solution S;

    cout << S.maximumProduct(nums);

    return 0;
}