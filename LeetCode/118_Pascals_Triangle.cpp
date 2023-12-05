#include <iostream>
#include <vector>
using namespace std;

class Solution 
{
public:
    vector<vector<int>> generate(int numRows) 
    {
        vector<vector<int>> vec;
        for(int i = 0; i < numRows; i++)
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

        return vec;
    }
};

void show(vector<vector<auto>> vec)
{
    for(int i = 0; i < vec.size(); i++)
    {
        for(int j = 0; j < vec[i].size(); j++)
        {
            cout << vec[i][j] << " ";
        }
        cout << endl;
    }
}

int main(int argc, char **argv)
{
    int numRows = 5;
    Solution S;
    vector<vector<int>> ans = S.generate(numRows);

    show(ans);

    return 0;
}