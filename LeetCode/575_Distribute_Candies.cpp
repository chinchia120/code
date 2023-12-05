#include <iostream>
#include <unordered_set>
#include <vector>
using namespace std;

class Solution 
{
public:
    int distributeCandies(vector<int>& candyType) 
    {
        unordered_set<int> dataset(candyType.begin(), candyType.end());
        //show_unordered_set(dataset);

        return min(dataset.size(), candyType.size()/2);
    }

    void show_unordered_set(unordered_set<int> set)
    {
        for(auto it: set) cout << it << " ";
        cout << '\n';
    }
};

int main(int argc, char ** argv)
{
    vector<int> candyType = {1, 1, 2, 2, 3, 3};
    Solution S;

    cout << S.distributeCandies(candyType);

    return 0;
}