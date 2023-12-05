#include <iostream>
#include <string>
#include <vector>
using namespace std;

class Solution 
{
public:
    vector<string> fizzBuzz(int n) 
    {   
        vector<string> vec;
        for(int i = 1; i <= n; i++)
        {   
            string str;
            if(i % 15 == 0)
            {
                str = "FizzBuzz";
            }
            else if(i % 5 == 0)
            {
                str = "Buzz";
            }
            else if(i % 3 == 0)
            {
                str = "Fizz";
            }
            else{
                str = to_string(i);
            }
            vec.push_back(str);
        }
        return vec;
    }

    void show_vector(vector<string> vec)
    {
        for(int i = 0; i < vec.size(); i++)
        {
            cout << vec[i] << " ";
        }
        cout << endl;
    }
};

int main(int argc, char **argv)
{
    int n = 15;
    Solution S;

    S.show_vector(S.fizzBuzz(n));

    return 0;
}