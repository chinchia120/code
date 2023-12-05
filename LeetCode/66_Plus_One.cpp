#include <iostream>
#include <vector>
using namespace std;

class Solution 
{
public:
    vector<int> plusOne(vector<int>& digits) 
    {   
        int flag = 0, tmp = 0, cnt = 0;
        vector<int> ans;
        
        for(int i = digits.size()-1; i >= 0; i--)
        {   
            if(flag == 0)
            {   
                flag = 1;
                tmp = digits[i] + 1;
            }
            else
            {
                tmp = digits[i] + cnt;
            }

            if(tmp >= 10)
            {
                cnt = 1;
                tmp = tmp - 10;
            }
            else
            {
                cnt = 0;
                tmp = tmp;
            }

            if(cnt == 1 && i == 0)
            {
                ans.insert(ans.begin(), tmp);
                ans.insert(ans.begin(), 1);
            }
            else
            {
                ans.insert(ans.begin(), tmp);
            }
            
        }

        return ans;
    }
};

int main(int argc, char **argv)
{
    vector<int> demo = {9   , 9, 9, 9};
    Solution S;
    vector<int> ans = S.plusOne(demo);
    
    for(int i = 0; i < ans.size(); i++)
    {
        cout << ans[i] << " ";
    }

    return 0;
}