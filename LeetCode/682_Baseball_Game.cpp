#include <iostream>
#include <vector>
#include <string>
#include <numeric>
using namespace std;

class Solution 
{
public:
    int calPoints(vector<string>& operations) 
    {
        vector<int> score;
        for(int i = 0; i < operations.size(); i++)
        {
            if(operations[i] == "+") score.push_back(score[score.size()-1]+score[score.size()-2]);
            else if(operations[i] == "D") score.push_back(score[score.size()-1]*2);
            else if(operations[i] == "C") score.erase(score.end()-1);
            else score.push_back(stoi(operations[i]));
            //show_vector(score);
        }
        
        return accumulate(score.begin(), score.end(), 0);
    }
    
    void show_vector(vector<int> vec)
    {
        for(int i = 0; i < vec.size(); i++) cout << vec[i] << " ";
        cout << endl;
    }
};

int main(int argc, char **argv)
{
    vector<string> operations = {"5","2","C","D","+"};
    Solution S;

    cout << S.calPoints(operations);

    return 0;
}