#include <iostream>
#include <vector>
using namespace std;

class Solution 
{
public:
    bool canPlaceFlowers(vector<int>& flowerbed, int n) 
    {   
        if(n == 0) return true;

        int sum = 0;
        for(int i = 0; i < flowerbed.size(); i++)
        {
            if(flowerbed[i] == 0)
            {   
                if(i-1 < 0 && i+1 >= flowerbed.size())
                {
                    sum++, i++;
                }
                else if(i-1 < 0)
                {
                    if(flowerbed[i+1] == 0) sum++, i++;
                    else continue;
                }
                else if(i+1 >= flowerbed.size())
                {
                    if(flowerbed[i-1] == 0) sum++, i++;
                    else continue;
                }
                else
                {
                    if(flowerbed[i-1] == 0 && flowerbed[i+1] == 0) sum++, i++;
                    else continue;
                }
                
            }
            if(n == sum) return true;
        }
        //cout << sum << endl;
        return false;
    }
};

int main(int argc, char **argv)
{
    vector<int> flowerbed = {0,1,0,1,0,1,0,0};
    int n = 1;
    Solution S;

    cout << S.canPlaceFlowers(flowerbed, n);

    return 0;
}