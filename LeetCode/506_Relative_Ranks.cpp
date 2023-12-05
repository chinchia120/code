#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
using namespace std;

class Solution 
{
public:
    vector<string> findRelativeRanks(vector<int>& score) 
    {
        vector<int> score2 = score;
        sort(score2.begin(), score2.end());
        reverse(score2.begin(), score2.end());
        //show_vector(score);
        //show_vector(score2);

        vector<string> rank;
        for(int i = 0; i < score.size(); i++)
        {
            for(int j = 0; j < score2.size(); j++)
            {
                if(score[i] == score2[j])
                {
                    if(j == 0) rank.push_back("Gold Medal");
                    else if(j == 1) rank.push_back("Silver Medal");
                    else if(j == 2) rank.push_back("Bronze Medal");
                    else rank.push_back(to_string(j+1));
                }
            }
        }

        return rank;
    }

    void show_vector(vector<int> vec)
    {
        for(int i = 0; i < vec.size(); i++) cout << vec[i] << " ";
        cout << endl;
    }
};

int main(int argc, char **argv)
{
    vector<int> score = {10, 3, 8, 9, 4};
    Solution S;

    S.show_vector(S.findRelativeRanks(score));

    return 0;
}