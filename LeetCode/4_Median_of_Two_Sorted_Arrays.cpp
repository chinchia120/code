#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

class Solution 
{
public:
    double findMedianSortedArrays(vector<int>& nums1, vector<int>& nums2) 
    {
        vector<int> list = nums1;
        for(int i = 0; i < nums2.size(); i++) list.push_back(nums2[i]);

        sort(list.begin(), list.end());

        if(list.size()%2 == 0) return (double)(list[list.size()/2]+list[list.size()/2-1])/2;
        else return (double)list[list.size()/2];
    }

    void show_vector(vector<int> vec)
    {
        for(int i = 0; i < vec.size(); i++) cout << vec[i] << " ";
        cout << endl;
    }
};

int main(int argc, char **argv)
{
    vector<int> nums1 = {1, 4};
    vector<int> nums2 = {2, 3};
    Solution S;

    cout << S.findMedianSortedArrays(nums1, nums2);

    return 0;
}