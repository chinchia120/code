#include <iostream>
#include <string>
using namespace std;
int main(void){
    string str;
    while(cin>>str){
        for(int i=0;i<str.length()+1;i++){
            int sum=1;
            for(int j=i+1;j<str.length()+1;j++){
                if(str[i]==str[j]){
                    sum+=1;
                }else{
                    if(sum>2){
                        cout<<sum<<str[i];
                    }else if(sum==2){
                        cout<<str[i]<<str[i];
                    }else{
                        cout<<str[i];
                    }
                    i=j-1;
                    break;
                }
            }
        }
        cout<<endl;
    }
}