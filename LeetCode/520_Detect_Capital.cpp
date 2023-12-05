#include <iostream>
#include <string>
using namespace std;

class Solution 
{
public:
    bool detectCapitalUse(string word) 
    {
        if(word.length() == 1) return true;

        int flag = 1, condition = 0;
        if(isupper(word[0])) condition = 1;
        if(islower(word[0])) condition = 2;

        for(int i = 1; i < word.length(); i++)
        {   
            if(condition == 1)
            {
                if(isupper(word[i])) condition = 3;
                if(islower(word[i])) condition = 4;
            }
            else if(condition == 2)
            {
                if(isupper(word[i])) flag = 0;
            }
            else if(condition == 3)
            {
                if(islower(word[i])) flag = 0;
            }
            else if(condition == 4)
            {
                if(isupper(word[i])) flag = 0;
            }

            if(flag == 0) return false;
        }    

        return true;
    }
};

int main(int argc, char **argv)
{
    string word = "USA";
    Solution S;

    cout << S.detectCapitalUse(word);

    return 0;
}