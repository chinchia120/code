#include <iostream>
#include <vector>
using namespace std;

class Solution 
{
public:
    int maxCount(int m, int n, vector<vector<int>>& ops) 
    {
        //show_2d_vector(ops);
        int min_m = m, min_n = n;
        for(int i = 0; i < ops.size(); i++)
        {
            if(ops[i][0] < min_m) min_m = ops[i][0];
            if(ops[i][1] < min_n) min_n = ops[i][1];
        }

        return min_m * min_n;
    }

    void show_2d_vector(vector<vector<int>> vec)
    {
        for(int i = 0; i < vec.size(); i++)
        {
            for(int j = 0; j < vec[i].size(); j++) cout << vec[i][j] << " ";
            cout << endl;
        }
    }
};

int main(int argc, char **argv)
{
    int m = 3, n = 3;
    vector<vector<int>> ops = {{2,2}, {3,3}, {3,3}, {3,3}, {2,2}, {3,3}, {3,3}, {3,3}, {2,2}, {3,3}, {3,3}, {3,3}};
    Solution S;

    cout << S.maxCount(m, n, ops);

    return 0;
}
