#include <iostream>
#include <vector>
#include <string>
using namespace std;

class Solution 
{
public:
    vector<int> selfDividingNumbers(int left, int right) 
    {
        vector<int> NumberList;
        for(int i = left; i <= right; i++)
        {
            if(isDividingNumbers(i)) NumberList.push_back(i);
        }
        return NumberList;
    }

    void show_vector(vector<int> vec)
    {
        for(int i = 0; i < vec.size(); i++) cout << vec[i] << " ";
        cout << endl;
    }

    bool isDividingNumbers(int num)
    {
        string str = to_string(num);
        for(int i = 0; i < str.length(); i++)
        {   
            if(str[i] == '0') return false;
            if(num%(str[i]-'0') != 0) return false;
        }
        return true;
    }
};

int main(int argc, char ** argv)
{
    int left = 1, right = 22;
    Solution S;

    S.show_vector(S.selfDividingNumbers(left, right));

    return 0;
}