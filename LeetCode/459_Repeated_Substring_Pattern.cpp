#include <iostream>
#include <string>
using namespace std;

class Solution 
{
public:
    bool repeatedSubstringPattern(string s) 
    {
        string str_1, str;
        str_1.push_back(s[0]);
        str.push_back(s[0]);
        for(int i = 1; i < s.length(); i++)
        {
            if(s[i] == str_1)
            {
                break;
            }
            else
            {
                str.push_back(s[i]);
            }
        }
        cout << str << endl;

        return true;
    }
};

int main(int argc, char **argv)
{
    string s = "abcabcabcabc";
    Solution S;
    cout << S.repeatedSubstringPattern(s);
}