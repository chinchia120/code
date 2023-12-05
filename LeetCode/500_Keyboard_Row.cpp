#include <iostream>
#include <string>
#include <vector>
#include <cctype>

using namespace std;

class Solution 
{
public:
    vector<string> findWords(vector<string>& words) 
    {
        //show_vector(words);
        vector<string> ans;

        for(int i = 0; i < words.size(); i++)
        {   
            
            string tmp = words[i];
            char tmp_0 = (char)tolower(tmp[0]);
            
            int flag = 0;
            if(tmp_0 == 'q' || tmp_0 == 'w' || tmp_0 == 'e' || tmp_0 == 'r' || tmp_0 == 't' || tmp_0 == 'y' || tmp_0 == 'u' || tmp_0 == 'i' || tmp_0 == 'o' || tmp_0 == 'p')
            {
                flag = 1;
            }
            else if(tmp_0 == 'a' || tmp_0 == 's' || tmp_0 == 'd' || tmp_0 == 'f' || tmp_0 == 'g' || tmp_0 == 'h' || tmp_0 == 'j' || tmp_0 == 'k' || tmp_0 == 'l')
            {
                flag = 2;
            }
            else
            {
                flag = 3;
            }

            int check = 1; 
            for(int i = 1; i < tmp.length(); i++)
            {   
                char tmp_i = (char)tolower(tmp[i]);
                if(flag == 1)
                {
                    if(tmp_i == 'q' || tmp_i == 'w' || tmp_i == 'e' || tmp_i == 'r' || tmp_i == 't' || tmp_i == 'y' || tmp_i == 'u' || tmp_i == 'i' || tmp_i == 'o' || tmp_i == 'p')
                    {
                        continue;
                    }
                    else
                    {
                        check = 0;
                        break;
                    }
                }
                
                if(flag == 2)
                {
                    if(tmp_i == 'a' || tmp_i == 's' || tmp_i == 'd' || tmp_i == 'f' || tmp_i == 'g' || tmp_i == 'h' || tmp_i == 'j' || tmp_i == 'k' || tmp_i == 'l')
                    {
                        continue;
                    }
                    else
                    {
                        check = 0;
                        break;
                    }
                }

                if(flag == 3)
                {
                    if(tmp_i == 'z' || tmp_i == 'x' || tmp_i == 'c' || tmp_i == 'v' || tmp_i == 'b' || tmp_i == 'n' || tmp_i == 'm')
                    {
                        continue;
                    }
                    else
                    {
                        check = 0;
                        break;
                    }
                }

                if(check == 0)
                {
                    break;
                }
            }

            if(check == 1)
            {
                ans.push_back(words[i]);
            }
        }

        return ans;
    }

    void show_vector(vector<string> vec)
    {
        for(int i = 0; i < vec.size(); i++)
        {
            cout << vec[i] << " ";
        }
        cout << endl;
    }
};

int main(int argc, char **argv)
{
    vector<string> words = {"Hello", "Alaska", "Dad", "Peace"};
    Solution S;
    
    S.show_vector(S.findWords(words));
    
    return 0;
}