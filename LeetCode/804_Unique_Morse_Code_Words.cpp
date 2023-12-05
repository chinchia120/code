#include <iostream>
#include <string>
#include <vector>
#include <unordered_set>
using namespace std;

class Solution
{
public:
    int uniqueMorseRepresentations(vector<string>& words)
    {
        unordered_set<string> morse_codes;
        for(int i = 0; i < words.size(); i++)
        {   
            string morse_code = "";
            for(int j = 0; j < words[i].size(); j++)
            {   
                if(words[i][j] == 'a') morse_code = push_back_char(morse_code, ".-");
                if(words[i][j] == 'b') morse_code = push_back_char(morse_code, "-...");
                if(words[i][j] == 'c') morse_code = push_back_char(morse_code, "-.-.");
                if(words[i][j] == 'd') morse_code = push_back_char(morse_code, "-..");
                if(words[i][j] == 'e') morse_code = push_back_char(morse_code, ".");
                if(words[i][j] == 'f') morse_code = push_back_char(morse_code, "..-.");
                if(words[i][j] == 'g') morse_code = push_back_char(morse_code, "--.");
                if(words[i][j] == 'h') morse_code = push_back_char(morse_code, "....");
                if(words[i][j] == 'i') morse_code = push_back_char(morse_code, "..");
                if(words[i][j] == 'j') morse_code = push_back_char(morse_code, ".---");
                if(words[i][j] == 'k') morse_code = push_back_char(morse_code, "-.-");
                if(words[i][j] == 'l') morse_code = push_back_char(morse_code, ".-..");
                if(words[i][j] == 'm') morse_code = push_back_char(morse_code, "--");
                if(words[i][j] == 'n') morse_code = push_back_char(morse_code, "-.");
                if(words[i][j] == 'o') morse_code = push_back_char(morse_code, "---");
                if(words[i][j] == 'p') morse_code = push_back_char(morse_code, ".--.");
                if(words[i][j] == 'q') morse_code = push_back_char(morse_code, "--.-");
                if(words[i][j] == 'r') morse_code = push_back_char(morse_code, ".-.");
                if(words[i][j] == 's') morse_code = push_back_char(morse_code, "...");
                if(words[i][j] == 't') morse_code = push_back_char(morse_code, "-");
                if(words[i][j] == 'u') morse_code = push_back_char(morse_code, "..-");
                if(words[i][j] == 'v') morse_code = push_back_char(morse_code, "...-");
                if(words[i][j] == 'w') morse_code = push_back_char(morse_code, ".--");
                if(words[i][j] == 'x') morse_code = push_back_char(morse_code, "-..-");
                if(words[i][j] == 'y') morse_code = push_back_char(morse_code, "-.--");
                if(words[i][j] == 'z') morse_code = push_back_char(morse_code, "--..");
            }
            morse_codes.insert(morse_code);
        }
        return morse_codes.size();
    }

    string push_back_char(string init_str, string add_str)
    {
        for(int i = 0; i < add_str.size(); i++) init_str.push_back(add_str[i]);
        return init_str;
    }
};

int main(int argc, char **argv)
{
    vector<string> words = {"gin","zen","gig","msg"};
    Solution S;

    cout << S.uniqueMorseRepresentations(words);

    return 0;
}