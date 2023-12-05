#include <iostream>
#include <string>
using namespace std;

class Solution 
{
public:
    string convertToBase7(int num) 
    {   
        string opt = "";
        if(num < 0){
            num = -num;
            opt = "-";
        }

        if(abs(num) < 7) return(opt+to_string(num));

        string str = "";
        while(num >= 7)
        {   
            
            str.insert(str.begin(), (num%7)+'0');
            num /= 7;

            if(num < 7)
            {   
                str.insert(str.begin(), num+'0');
            }
        }
        return opt+str;
    }
};

int main(int argc, char **argv)
{
    int num = 0;
    Solution S;

    cout << S.convertToBase7(num);

    return 0;
}