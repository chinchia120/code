#include <iostream>
#include <string>
#include <vector>
using namespace std;

class Solution 
{
public:
    bool isValid(string s) 
    {   
        if(s.length()%2 != 0)
        {
            return false;
        }
        vector<char> str;
        for(int i = 0; i < s.length(); i++)
        {
            str.push_back(s.at(i));
        }

        int flag = 0;
        int str_len = str.size();
        while(str_len != 0)
        {
            int index = str_len/2 - 1;
            //cout << str[index] << " " << str[index+1] << endl;

            for(int i = 0; i < str.size(); i++)
            {
                cout << str[i] << " ";
            }
            cout << endl;

            if(str[index] == '(' && str[index+1] == ')')
            {
                auto del = str.erase(str.begin()+str.size()/2-1, str.begin()+str.size()/2+1);

                str_len = str.size();
            }
            else if(str[index] == '{' && str[index+1] == '}')
            {
                auto del = str.erase(str.begin()+str.size()/2-1, str.begin()+str.size()/2+1);

                str_len = str.size();
            }
            else if(str[index] == '[' && str[index+1] == ']')
            {
                auto del = str.erase(str.begin()+str.size()/2-1, str.begin()+str.size()/2+1);

                str_len = str.size();
            }
            else
            {
                flag = 1;
                break;
            }
        }

        if(flag == 0)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
};

int main(int argc, char **argv)
{
    string str = "()[]{}";
    Solution S;

    bool ans = S.isValid(str);
    cout << ans;
    return 0;
}