#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

class Solution 
{
public:
    vector<int> intersect(vector<int>& nums1, vector<int>& nums2) 
    {
        sort(nums1.begin(), nums1.end());
        sort(nums2.begin(), nums2.end());

        vector<int> ans;
        int index = 0;
        for(int i = 0; i < nums1.size(); i++)
        {
            for(int j = index; j < nums2.size(); j++)
            {
                if(nums1[i] == nums2[j])
                {   
                    ans.push_back(nums1[i]);
                    index = j+1;
                    break;
                }
            }
        }

        return ans;
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
    vector<int> nums1 = {4,3,9,3,1,9,7,6,9,7};
    vector<int> nums2 = {5,0,8};
    Solution S;

    vector<int> ans = S.intersect(nums1, nums2);
    S.show_vector(ans);

    return 0;
}