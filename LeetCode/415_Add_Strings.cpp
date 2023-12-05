#include <iostream>
#include <string>
using namespace std;

class Solution 
{
public:
    string addStrings(string num1, string num2) 
    {   
        if(num1.length() < num2.length())
        {
            string tmp = num1;
            num1 = num2;
            num2 = tmp;
        }
        
        num1 = inverse_string(num1);
        num2 = inverse_string(num2);
        for(int i = 0; i < num1.length(); i++)
        {
            if(i > num2.length()-1)
            {
                num2.push_back('0');
            }        
        }
        num1 = inverse_string(num1);
        num2 = inverse_string(num2);

        string sum = "";
        int next = 0;
        for(int i = 0; i < num1.length(); i++)
        {   
            int dig_num1 = num1.at(num1.length()-1-i)-'0', dig_num2 = num2.at(num2.length()-1-i)-'0';
            if(i == num1.length()-1)
            {   
                sum.insert(0, to_string(dig_num1+dig_num2+next));
                break;
            }

            if(dig_num1+dig_num2+next >= 10)
            {   
                sum.insert(0, to_string(dig_num1+dig_num2+next-10));
                next = 1;
            }
            else
            {
                sum.insert(0, to_string(dig_num1+dig_num2+next));
                next = 0;
            }
        }
       return sum;
    }

    string inverse_string(string str)
    {   
        string inv_str = "";
        for(int i = 0; i < str.length(); i++)
        {
            inv_str.push_back(str.at(str.length()-1-i));
        }
        return inv_str;
    }
};

int main(int argc, char **argv)
{
    //string num1 = "123456789" , num2 = "987654321";
    //string num1 = "11" , num2 = "123";
    //string num1 = "10" , num2 = "90";
    string num1 = "584", num2 = "18";
    Solution S;

    cout << S.addStrings(num1, num2);

    return 0;
}