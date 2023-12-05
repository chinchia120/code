#include <iostream>
#include <vector>
using namespace std;

class Solution 
{
public:
    int islandPerimeter(vector<vector<int>>& grid) 
    {   
        //show_2d_vector(grid);
        
        int sum = 0, edge_i = grid.size()-1, edge_j = grid[0].size()-1;
        for(int i = 0; i < grid.size(); i++)
        {
            for(int j = 0; j < grid[i].size(); j++)
            {
                if(grid[i][j] == 1)
                {   
                    //cout << "index = " << i << ", " << j << endl;
                    if(check_index(i, j, -1, 0, edge_i, edge_j)) 
                    {
                        if(grid[i-1][j+0] == 0) sum++;
                    }
                    else
                    {
                        sum++;
                    }

                    if(check_index(i, j, +1, 0, edge_i, edge_j)) 
                    {
                        if(grid[i+1][j+0] == 0) sum++;
                    } 
                    else
                    {
                        sum++;
                    }

                    if(check_index(i, j, 0, -1, edge_i, edge_j)) 
                    {
                        if(grid[i+0][j-1] == 0) sum++;
                    } 
                    else
                    {
                        sum++;
                    }

                    if(check_index(i, j, 0, +1, edge_i, edge_j)) 
                    {
                        if(grid[i+0][j+1] == 0) sum++;
                    }
                    else
                    {
                        sum++;
                    }
                    
                }
            }
        }

        return sum;
    }

    bool check_index(int i, int j, int opt_i, int opt_j, int max_i, int max_j)
    {   
        //cout << i + opt_i << " " << j + opt_j << endl;
        if(i + opt_i < 0 || i + opt_i > max_i || j + opt_j < 0 || j + opt_j > max_j) return false;
        else return true; 
    }

    void show_2d_vector(vector<vector<int>> vec)
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
};

int main(int argc, char **argv)
{   
    vector<vector<int>> grid = {{0,1,0,0}, {1,1,1,0}, {0,1,0,0}, {1,1,0,0}};
    Solution S;

    cout << S.islandPerimeter(grid);

    return 0;
}