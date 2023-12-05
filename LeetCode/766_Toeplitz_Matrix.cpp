#include <iostream>
#include <vector>
using namespace std;

class Solution
{
public:
    bool isToeplitzMatrix(vector<vector<int>>& matrix)
    {
        //show_2d_vector(matrix);

        for(int i = 1; i < matrix.size(); i++)
        {
            for(int j = 1; j < matrix[i].size(); j++)
            {
                if(matrix[i][j] != matrix[i-1][j-1]) return false;
            }
        }

        return true;
    }

    void show_2d_vector(vector<vector<int>>& vec)
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
    vector<vector<int>> matrix = {{1,2,3,4},{5,1,2,3},{9,5,1,2}};
    Solution S;

    cout << S.isToeplitzMatrix(matrix);

    return 0;
}