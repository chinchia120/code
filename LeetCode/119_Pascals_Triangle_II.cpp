#include <iostream>
#include <vector>
using namespace std;

class Solution 
{
public:
    vector<int> getRow(int rowIndex) 
    {
        vector<vector<int>> vec;
        for(int i = 0; i < rowIndex+1; i++)
        {   
            vec.push_back(vector<int>());
            for(int j = 0; j <= i; j++)
            {
                if(j == 0 || j == i)
                {
                    vec[i].push_back(1);
                }
                else
                {
                    vec[i].push_back(vec[i-1][j-1] + vec[i-1][j]);
                }   
            }
        }

        return vec[rowIndex];
    }
};

void show(vector<auto> vec)
{
    for(int i = 0; i < vec.size(); i++)
    {
        cout << vec[i] << " ";
    }
}

int main(int argc, char **argv)
{
    int rowIndex = 3;
    Solution S;
    vector<int> ans = S.getRow(rowIndex);

    show(ans);

    return 0;
}