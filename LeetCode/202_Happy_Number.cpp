#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

class Solution 
{
private:
    vector<int> check_number(int n)
    {   
        vector<int> my_num;
        while(true)
        {   
            if(n < 10)
            {
                my_num.push_back(n);

                break;
            }

            if(n > 9 && n < 100)
            {
                my_num.insert(my_num.begin(), n%10);
                my_num.insert(my_num.begin(), n/10);

                break;
            }
            else
            {
                my_num.insert(my_num.begin(), n%10);
                n /= 10;
            }
        }
        //show_vector(my_num);

        return my_num;
    }

    void show_vector(vector<auto> vec)
    {
        for(int i = 0; i < vec.size(); i++)
        {
            cout << vec[i] << " ";
        }
        cout << endl;
    }

public:
    bool isHappy(int n) 
    {   
        vector<int> num = check_number(n);
        vector<int> same = {n};

        while(true)
        {   
            int sum = 0;
            for(int i = 0; i < num.size(); i++)
            {
                sum += num[i]*num[i];
            }
            
            if(sum == 1)
            {
                return true;
            }
            else
            {   
                for(int j = 0; j < same.size(); j++)
                {
                    if(same[j] == sum)
                    {
                        return false;
                    }
                }
                same.push_back(sum);
            }

            num = check_number(sum);
        }

        return false;
    }
};

int main(int argc, char **argv)
{
    int n = 2;
    Solution S;
    bool ans = S.isHappy(n);

    cout << ans;

    return 0;
}