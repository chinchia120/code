#include <iostream>
#include <vector>
using namespace std;

class NumArray 
{
public:
    vector<int> vec;
    NumArray(vector<int>& nums) 
    {
        vec = nums;
    }
    
    int sumRange(int left, int right) 
    {   
        int sum = 0;
        for(int i = left; i <= right; i++)
        {
            sum += vec[i];
        }

        return sum;
    }
};

int main(int argc, char **argv)
{
    
    return 0;
}