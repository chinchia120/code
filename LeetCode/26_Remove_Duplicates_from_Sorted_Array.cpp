#include <iostream>
#include <vector>
using namespace std;

class Solution 
{
public:
    int removeDuplicates(vector<int>& nums) 
    {   
        int cnt = 0;
        int arr[] = {nums[0]};

        for(int i = 1; i < nums.size(); i++)
        {
            if(arr[cnt] != nums[i])
            {
                cnt++;
                arr[cnt] = nums[i];
            }
        }

        return arr[1];
    }
};

int main(int argc, char **argv)
{
    vector<int> demo = {0, 0, 1, 1, 1, 2, 2, 3, 3, 4};
    Solution S;
    int ans = S.removeDuplicates(demo);

    cout << ans << " ";

    return 0;
}