#include <iostream>
#include <string>
#include <vector>
using namespace std;

class Solution 
{
public:
    bool rotateString(string s, string goal) 
    {   
        if(s == goal) return true;

        vector<int> index = {};
        for(int i = 0; i < goal.length(); i++) if(s[0] == goal[i]) index.push_back(i);

        for(int i = 0; i < index.size(); i++)
        {   
            string check = "";
            for(int j = 0; j < goal.size(); j++)
            {
                if(index[i]+j < goal.size()) check.push_back(goal[index[i]+j]);
                else check.push_back(goal[index[i]+j-goal.size()]);
            }
            //cout << s << " " << check << endl;
            if(check == s) return true;
        }
        return false;
    }
};

int main(int argc, char **argv)
{
    string s = "abcde", goal = "cdeab";
    Solution S;

    cout << S.rotateString(s, goal);

    return 0;
}