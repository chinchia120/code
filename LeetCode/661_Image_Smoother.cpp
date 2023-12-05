#include <iostream>
#include <vector>
using namespace std;

class Solution
{
public:
    vector<vector<int>> imageSmoother(vector<vector<int>>& img)
    {   
        vector<vector<int>> img_smooth = img;
        for(int i = 0; i < img.size(); i++)
        {
            for(int j = 0; j < img[i].size(); j++)
            {
                int sum = 0, cnt = 0;
                for(int row = -1; row < 2; row++)
                {
                    for(int column = -1; column < 2; column++)
                    {
                        if(i+row >= 0 && j+column >= 0 && i+row < img.size() && j+column < img[i].size())
                        {   
                            sum += img[i+row][j+column];
                            cnt++;
                        } 
                    }
                }
                img_smooth[i][j] = sum/cnt;
            }
        }

        return img_smooth;
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
    vector<vector<int>> img = {{100,200,100}, {200,50,200}, {100,200,100}};
    Solution S;

    S.show_2d_vector(S.imageSmoother(img));

    return 0;
}