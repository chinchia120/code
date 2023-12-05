#include <iostream>
#include <string>
#include <algorithm>
using namespace std;

class Solution 
{
public:
    char findTheDifference(string s, string t) 
    {   
        sort(s.begin(), s.end());
        sort(t.begin(), t.end());

        for(int i = 0; i < s.length(); i++)
        {
            if(s.at(i) != t.at(i))
            {   
                return t.at(i);
            }
        }
        
        return t.at(t.length()-1);
    }
};


int main(int argc, char **argv)
{
    string s = "abcd";
    string t = "abcde";
    Solution S;
    
    cout << S.findTheDifference(s, t);

    return 0;
}