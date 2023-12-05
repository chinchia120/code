#include <iostream>
#include <string>
using namespace std;

class Solution 
{
public:
    string reverseVowels(string s) 
    {   
        int index = s.length()-1;
        int i = 0, j = index;
        for(i = 0; i < s.length(); i++)
        {   
            for(j = index; j > i; j--)
            {   
                if((tolower(s.at(i)) == 'a' || tolower(s.at(i)) == 'e' || tolower(s.at(i)) == 'i' || tolower(s.at(i)) == 'o' || tolower(s.at(i)) == 'u') && (tolower(s.at(j)) == 'a' || tolower(s.at(j)) == 'e' || tolower(s.at(j)) == 'i' || tolower(s.at(j)) == 'o' || tolower(s.at(j)) == 'u'))
                {   
                    char tmp = s.at(j);
                    s.at(j) = s.at(i);
                    s.at(i) = tmp;
                    index = j-1;
                    break;
                }
            }

            if(i > j)
            {
                break;
            }
        }
        return s;
    }
};

int main(int argc, char **argv)
{
    string s = "leetcode";
    Solution S;

    cout << S.reverseVowels(s) << endl;

    return 0;
}
