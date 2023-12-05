#include <iostream>
#include <string>
#include <vector>
using namespace std;

class Solution 
{
public:
    void reverseString(vector<char>& s) 
    {
        for(int i = 0; i < s.size()/2; i++)
        {
            swap(s[i], s[s.size()-1-i]);
        }
    }
};

int main(int argc, char **argv)
{
    vector<char> s = {'h', 'e', 'l', 'l', 'o'};
    Solution S;

    S.reverseString(s);

    return 0;
}