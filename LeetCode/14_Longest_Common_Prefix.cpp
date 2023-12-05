#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;

class Solution 
{
public:
    string longestCommonPrefix(vector<string>& strs) 
    {
        string LongestString = "";

        for(int i = 0; i < strs[0].length(); i++)
        {
            string str = strs[0];
            char ch = str.at(i);

            int flag = 1;
            for(int j = 1; j < strs.size(); j++)
            {   
                if((i+1) > strs[j].length())
                {
                    flag = 0;
                    break;
                }

                string str_check = strs[j];
                char ch_check = str_check.at(i);

                if(ch != ch_check)
                {
                    flag = 0;
                    break;
                }
            }

            if(flag == 1) LongestString.push_back(ch);
            if(flag == 0) break;
        }

        return LongestString;
    }

    void show_vector(vector<string> vec)
    {
        for(int i = 0; i < vec.size(); i++) cout << vec[i] << " ";
        cout << "\n";
    }
};

int main(int argc, char **argv)
{
    vector<string> strs = {"flower", "flow", "flowight"};
    Solution S;

    cout << S.longestCommonPrefix(strs);

    return 0;
}