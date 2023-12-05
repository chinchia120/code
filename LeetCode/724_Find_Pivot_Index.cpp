#include <iostream>
#include <vector>
using namespace std;

class Solution 
{
public:
    int pivotIndex(vector<int>& nums) 
    {
        for(int i = 0; i < nums.size(); i++)
        {
            int sumR = 0, sumL = 0;
            for(int j = 0; j < i; j++) sumL += nums[j];
            for(int j = nums.size()-1; j > i; j--) sumR += nums[j];
            if(sumR == sumL) return i;
        }
        return -1;
    }
};

int main(int argc, char **argv)
{
    vector<int> nums = {1,7,3,6,5,6};
    Solution S;

    cout << S.pivotIndex(nums);

    return 0;
}