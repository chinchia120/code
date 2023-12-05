#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

class Solution
{
public:
    double largestTriangleArea(vector<vector<int>>& points)
    {
        //show_2d_vector(points);
        double area = 0;
        for(int i = 0; i < points.size(); i++)
        {
            for(int j = i+1; j < points.size(); j++)
            {
                for(int k = j+1; k < points.size(); k++)
                {
                    int x1 = points[j][0] - points[i][0];
                    int x2 = points[k][0] - points[i][0];
                    int y1 = points[j][1] - points[i][1];
                    int y2 = points[k][1] - points[i][1];
                    if(abs(x1*y2-x2*y1)*0.5 > area) area = abs(x1*y2-x2*y1)*0.5;
                }
            }
        }
        return area;
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
    vector<vector<int>> points = {{0,0},{0,1},{1,0},{0,2},{2,0}};
    Solution S;

    cout << S.largestTriangleArea(points);

    return 0;
}