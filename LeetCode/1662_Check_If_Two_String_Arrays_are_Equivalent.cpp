#include <iostream>
#include <vector>
using namespace std;

class Solution
{
public:
    bool arrayStringsAreEqual(vector<string>& word1, vector<string>& word2)
    {
        string str1 = "", str2 = "";
        for(int i = 0; i < word1.size(); i++)
        {
            for(int j = 0; j < word1[i].size(); j++) str1.push_back(word1[i][j]);
        }

        for(int i = 0; i < word2.size(); i++)
        {
            for(int j = 0; j < word2[i].size(); j++) str2.push_back(word2[i][j]);
        }

        if(str1 == str2) return true;
        else return false;
    }
};

int main(int argc, char **argv)
{
    vector<string> word1 = {"abc", "d", "defg"};
    vector<string> word2 = {"abcddef"};
    Solution S;

    cout << S.arrayStringsAreEqual(word1, word2);

    return 0;
}