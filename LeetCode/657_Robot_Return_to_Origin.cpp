#include <iostream>
#include <string>
using namespace std;

class Solution 
{
public:
    bool judgeCircle(string moves) 
    {   
        int x = 0, y = 0;
        for(int i = 0; i < moves.length(); i++)
        {
            if(moves[i] == 'U') y++;
            if(moves[i] == 'D') y--;
            if(moves[i] == 'R') x++;
            if(moves[i] == 'L') x--;
            //cout << x << " " << y << endl;
        }
        
        if(x == 0 && y == 0) return true;
        else return false;
    }
};

int main(int argc, char **argv)
{
    string move = "UD";
    Solution S;

    cout << S.judgeCircle(move);
}