#include <iostream>
#include <string>
using namespace std;

class Solution 
{
public:
    string addBinary(string a, string b) 
    {   
        int ans = stoi(a+b);
        
        return to_string(ans);
    }
};

int main(int argc, char **argv)
{
    string a = "11", b = "1";
    Solution S;
    string ans = S.addBinary(a, b);

    cout << ans;

    return 0;
}