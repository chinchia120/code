#include <iostream>
#include <vector>
using namespace std;

int main(int argc, char **argv)
{
    vector<int> v1 = {1, 2, 3};

    v1.push_back(4);
    v1.insert(v1.begin(), 0);
    
    for(int i = 0; i < v1.size(); i++)
    {
        cout << v1[i] << " ";
    }
    cout << endl;

    for(auto i = v1.begin(); i != v1.end(); i++)
    {  
        cout << *i << " "; 
    }
    cout << endl;
    
    for(auto &i : v1)
    {
        cout << i << " ";
    }
    cout << endl;

    return 0;
}