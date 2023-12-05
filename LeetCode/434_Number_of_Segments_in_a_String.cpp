#include <iostream>
#include <string>
using namespace std;

class Solution 
{
public:
    int countSegments(string s) 
    {   
        if(s == "")
        {
            return 0;
        }

        int cnt = 0;
        for(int i = 0; i < s.length()-1; i++)
        {
            if((int)s.at(i) != 32 && (int)s.at(i+1) == 32) cnt++;
        }

        if((int)s.at(s.length()-1) != 32) cnt++;
        return cnt;
    }
};

int main(int argc, char **argv)
{
    string s = "Of all the gin joints in all the towns in all the world,   ";
    Solution S;

    cout << S.countSegments(s);

    return 0;
}