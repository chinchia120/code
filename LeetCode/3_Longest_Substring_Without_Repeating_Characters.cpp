#include <iostream>
#include <string>
#include <vector>
using namespace std;

class Solution 
{
public:
    int lengthOfLongestSubstring(string s) 
    {
        vector<char> cnt_str;
        cnt_str.push_back(s.at(0));

        for(int i = 1; i < s.length(); i++)
        {   
            int flag = 0;
            for(int j = 0; j < cnt_str.size(); j++)
            {   
                if(s.at(i) == cnt_str[j])
                {   
                    flag = 1;
                    break;
                }        
            }
            if(flag == 0)
            {
                cnt_str.push_back(s.at(i));
            }
        }
        return cnt_str.size();
    }
};

int main(int argc, char **argv)
{
    string str = "pwwkew";

    Solution S;
    int ans = S.lengthOfLongestSubstring(str);

    cout << ans;

    return 0;
}