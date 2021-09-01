#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <stack>
#include <algorithm>
#include <vector>
using namespace std;

int main() {
    int n1;
    cin>>n1;
    vector<string> name[n1];
    for(int i=0;i<n1;i++){
        string _name;
        cin>>_name;
        name[i].push_back(_name);
    }
    int n2;
    cin>>n2;
    for(int i=0;i<n2;i++){
        string _name1,_name2;
        cin>>_name1>>_name2;
        for(int j=0;j<n1;j++){
            if(name[j][0]==_name1){
                name[j].push_back(_name2);
            }
            if(name[j][0]==_name2){
                name[j].push_back(_name1);
            }
        }
    }
    for(int i=0;i<n1;i++){
        for(int j=1;j<name[i].size();j++){
            for(int k=j+1;k<name[i].size();k++){
                if(name[i][j]>name[i][k]){
                    string tmp=name[i][j];
                    name[i][j]=name[i][k];
                    name[i][k]=tmp;
                }
            }
        }
    }
    for(int i=0;i<n1;i++){
        for(int j=0;j<name[i].size();j++){
            if(j==0){
                cout<<name[i][j]<<":";
            }else{
                cout<<" "<<name[i][j];
            }
        }
        cout<<endl;
    }
    return 0;
}
