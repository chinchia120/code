#include <iostream>
#include <algorithm>
#include <string>
using namespace std;

class Solution 
{
public:
    bool isAnagram(string s, string t) 
    {
        sort(s.begin(), s.end());
        sort(t.begin(), t.end());

        if(s.compare(t) == 0) return true;
        else return false;
    }
};

int main(int argc, char **argv)
{
    string s = "anagram", t = "nagaram";
    Solution S;

    cout << S.isAnagram(s, t);

    return 0;
}