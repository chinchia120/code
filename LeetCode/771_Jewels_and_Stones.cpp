#include <iostream>
#include <string>
#include <algorithm>
using namespace std;

class Solution 
{
public:
    int numJewelsInStones(string jewels, string stones) 
    {
        sort(jewels.begin(), jewels.end());
        sort(stones.begin(), stones.end());

        int sum = 0, cnt = 0;
        for(int i = 0; i < jewels.length(); i++)
        {
            int jewel = (int)jewels[i];
            for(int j = cnt; j < stones.length(); j++)
            {
                int stone = (int)stones[j];
                //cout << jewel << " " << stone << endl;
                if(jewel == stone)
                {
                    sum++;
                }
                if(jewel < stone)
                {
                    cnt = j;
                    break;
                }
            }
        }
        return sum;
    }
};

int main(int atgc, char **argv)
{
    string jewels = "dbca", stones = "beb";
    Solution S;

    cout << S.numJewelsInStones(jewels, stones);

    return 0;
}