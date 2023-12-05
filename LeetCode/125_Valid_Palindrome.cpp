#include <iostream>
#include <string>
#include <algorithm>
using namespace std;

class Solution 
{
public:
    bool isPalindrome(string s) 
    {   
        for(int i = 0; i < s.length(); i++)
        {
            s[i] = tolower(s[i]);
            if((int)s[i] < 48 || (int)s[i] > 122 || ((int)s[i] > 57 && (int)s[i] < 97))
            {
                s.replace(i, 1, " ");
            }
        }
        s.erase(remove(s.begin(), s.end(), ' '), s.end());

        int flag = 0;
        for(int i = 0; i < (int)s.length()/2; i++)
        {
            if(s[i] != s[s.length()-1-i])
            {
                flag = 1;
            }
        }
        
        if(flag == 0)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
};

int main(int argc, char **argv)
{
    string s = "0P";
    Solution S;
    bool ans = S.isPalindrome(s);

    cout << ans;

    return 0;
}