#include <iostream>
#include <vector>
using namespace std;

class Solution 
{
public:
    void moveZeroes(vector<int>& nums) 
    {
        vector<int> index;
        for(int i = 0; i < nums.size(); i++)
        {
            for(int j = i+1; j < nums.size(); j++)
            {
                if(nums[i] == 0)
                {
                    int tmp = nums[j];
                    nums[j] = nums[i];
                    nums[i] = tmp;
                }
            }
        }
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
    vector<int> nums = {-1, 0, -2 , 0, 1, 3, 0};
    Solution S;

    S.moveZeroes(nums);
    S.show_vector(nums);

    return 0;
}