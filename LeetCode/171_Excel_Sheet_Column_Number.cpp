#include <iostream>
#include <string>
#include <vector>
#include <cmath>
using namespace std;

class Solution 
{
public:
    int titleToNumber(string columnTitle) 
    {   
        vector<int> vec;
        int sum = 0;

        for(int i = 0; i < columnTitle.length(); i++)
        {
            vec.push_back(((int)columnTitle[i]+'0')-112);
        }

        for(int i = 0; i < vec.size(); i++)
        {
            sum += pow(26, vec.size()-1-i)*vec[i];
        }

        return sum;
    }
};

int main(int argc, char **argv)
{
    string columnTitle = "AA";
    Solution S;
    int ans = S.titleToNumber(columnTitle);

    cout << ans;

    return 0;
}