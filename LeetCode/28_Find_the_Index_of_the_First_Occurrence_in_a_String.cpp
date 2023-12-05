#include <iostream>
#include <string>
using namespace std;

class Solution 
{
public:
    int strStr(string haystack, string needle) 
    {   
        int i = 0, j = 0, flag = 0;
        
        for(i = 0; i < haystack.length(); i++)
        {   
            flag = 0;

            for(j = 0; j < needle.length(); j++)
            {
                if(haystack[i+j] != needle[j])
                {
                    flag = 1;
                    break;
                }
            }
            
            if(flag == 0)
            {
                return i;
            }
        }

        return -1;
    }
};

int main(int argc, char **agrv)
{
    string haystack = "sadbutsad", needle = "sad";

    Solution S;
    int ans = S.strStr(haystack, needle);

    cout << ans;

    return 0;
}