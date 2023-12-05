#include <iostream>
#include <vector>
using namespace std;

class Solution 
{
public:
    int addDigits(int num) 
    {   
        while(num > 9)
        {   
            vector<int> vec;
            int tmp = num, sum = 0;
            while(tmp)
            {   
                vec.push_back(tmp%10);
                tmp /= 10;
            }

            for(int i = 0; i < vec.size(); i++)
            {
                sum += vec[i];
            }
            num = sum;
        }
        return num;
    }
};

int main(int argc, char **argv)
{
    int num = 10;
    Solution S;

    cout << S.addDigits(num);

    return 0;
}