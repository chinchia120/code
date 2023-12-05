#include <iostream>
#include <string>
using namespace std;

class Solution 
{
public:
    string toLowerCase(string s) 
    {
        string str;
        for(int i = 0; i < s.length(); i++) str.push_back((char)tolower(s[i]));
        return str;
    }
};

int main(int argc, char **argv)
{
    string s = "Hello";
    Solution S;

    cout << S.toLowerCase(s);

    return 0;
}