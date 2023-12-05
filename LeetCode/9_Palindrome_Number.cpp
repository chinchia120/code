#include <iostream>
#include <string>
using namespace std;

class Solution 
{
public:
    bool isPalindrome(int x) 
    {
        if(x < 0)
        {
            return false;
        }
        else
        {
            string str1 = to_string(x);
            
            int flag = 0;
            for(int i = 0; i < str1.length(); i++)
            {
                if(str1.at(i) != str1.at(str1.length()-1-i))
                {
                    flag = 1;
                    break;
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
    }
};

int main(int argc, char **argv)
{
    int x = 121;
    Solution S;

    bool ans = S.isPalindrome(x);
    cout << ans;
    return 0;
}