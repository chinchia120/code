#include <iostream>
#include <string>
#include <vector>
using namespace std;

class Solution 
{
public:
    string convertToTitle(int columnNumber) 
    {   
        vector<int> vec;
        int integer = columnNumber/26;
        int decimal = columnNumber%26;

        while(1)
        {    
            //cout << integer << " " << decimal << endl;

            if(decimal == 0)
            {
                integer -= 1;
                decimal = 26;
            }

            vec.push_back(decimal);

            if(integer < 27)
            {   
                vec.push_back(integer);
                break;
            }
            
            decimal = integer%26;
            integer = integer/26;
        }

        string str;
        for(int i = vec.size()-1; i >= 0; i--)
        {   
            //cout << vec[i] << " ";
            if(vec[i] != 0)
            {
                str.push_back(vec[i]+'@');
            }
        }
        
        return str;
    }
};

int main(int argc, char **argv)
{
    int columnNumber = 1000;
    Solution S;
    string ans = S.convertToTitle(columnNumber);

    cout << ans;

    return 0;
}