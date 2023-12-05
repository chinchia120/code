#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

class Solution
{
public:
    double findMaxAverage(vector<int>& nums, int k)
    {
        double avg, sum = 0;
        for(int i = 0; i < k; i++) sum += nums[i];
        avg = sum/k;
        for(int i = 1; i < nums.size()-k+1; i++)
        {
            sum = sum - nums[i-1] + nums[i+k-1];
            if(sum/k > avg) avg = sum/k;
        }

        return avg;
    }
};

int main(int argc, char ** argv)
{
    vector<int> nums = {1,12,-5,-6,50,3};
    int k = 4;
    Solution S;

    cout << S.findMaxAverage(nums, k);

    return 0;
}