#include <iostream>
#include <string>
using namespace std;

class Solution 
{
public:
    int lengthOfLastWord(string s) 
    {   
        int cnt_last_wold = 0, flag = 0;
        for(int i = s.length()-1; i >= 0; i--)
        {
            if(flag == 0 && s.at(i) != ' ')
            {
                flag = 1;
                cnt_last_wold = i;
            }
            
            if(flag == 1 && s.at(i) == ' ') 
            {
                return cnt_last_wold-i;
            }
        }
        
        if(flag == 1 && s.at(0) != ' ')
        {
            return cnt_last_wold+1;
        } 

        return s.length();
    }
};

int main(int argc, char **argv)
{
    string demo = "abc abc   ";
    Solution S;
    int ans = S.lengthOfLastWord(demo);

    cout << ans;

    return 0;
}