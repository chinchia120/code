#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

class Solution 
{
public:
    vector<int> intersection(vector<int>& nums1, vector<int>& nums2) 
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
        ans = remove_same(ans);

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

    vector<int> remove_same(vector<int> vec)
    {   
        if(vec.size() == 0)
        {
            return vec;
        }

        vector<int> rm_vec = {vec[0]};
        int cnt = 0;
        for(int i = 1; i < vec.size(); i++)
        {
            if(vec[i] != rm_vec[cnt])
            {
                rm_vec.push_back(vec[i]);
                cnt++;
            }
        }
        //show_vector(rm_vec);

        return rm_vec;
    }
};

int main(int argc, char **argv)
{
    vector<int> nums1 = {4,3,9,3,1,9,7,6,9,7};
    vector<int> nums2 = {5,0,8};
    Solution S;

    vector<int> ans = S.intersection(nums1, nums2);
    S.show_vector(ans);

    return 0;
}