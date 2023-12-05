#include <iostream>
#include <vector>
using namespace std;

class Solution
{
public:
    char nextGreatestLetter(vector<char>& letters, char target)
    {
        for(int i = 0; i < letters.size(); i++)
        {
            //cout << (int)letters[i] << " ";
            if((int)target < (int)letters[i]) return letters[i];
        }    

        return letters[0];
    }
};

int main(int argc, char **argv)
{
    vector<char> letters = {'c', 'f', 'j'};
    char target = 'a';
    Solution S;

    cout << S.nextGreatestLetter(letters, target);

    return 0;
}