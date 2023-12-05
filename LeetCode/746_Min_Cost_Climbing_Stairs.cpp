#include <iostream>
#include <vector>
using namespace std;

class Solution
{
public:
    int minCostClimbingStairs(vector<int>& cost)
    {   
        //show_vector(cost);
        for(int i = 0; i < 2; i++) cost.push_back(0);
        vector<int> sum_cost = {cost[0], cost[1]};
        for(int i = 2; i < cost.size(); i++)
        {   
            cout << sum_cost[i-1]+cost[i] << " " << sum_cost[i-2]+cost[i] << endl;
            sum_cost.push_back((min(sum_cost[i-1]+cost[i], sum_cost[i-2]+cost[i])));
        }
        //show_vector(sum_cost);
        
        return sum_cost.back();
    }

    void show_vector(vector<int> vec)
    {
        for(int i = 0; i < vec.size(); i++) cout << vec[i] << " ";
        cout << endl;
    }
};

int main(int argc, char **argv)
{
    //vector<int> cost = {1,100,1,1,1,100,1,1,100,1};
    vector<int> cost = {10, 15, 20};
    Solution S;

    cout << S.minCostClimbingStairs(cost);

    return 0;
}