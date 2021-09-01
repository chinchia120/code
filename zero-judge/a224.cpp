#include <iostream>
#include <string>
#include <cstdlib>
using namespace std;
int main(void){
    string str;
    while(cin>>str){
        string arr[str.length()];
        int k=0;
        for(int i=0;i<str.length();i++){
            if(isupper(str.at(i))){
                arr[k]=char(str.at(i)+32);
                //cout<<arr[k]<<" ";
                k++;
            }else if(islower(str.at(i))){
                arr[k]=str.at(i);
                //cout<<arr[k]<<" ";
                k++;
            }else{
                continue;
            }
        }
        /*for(int i=0;i<k;i++){
            cout<<arr[i]<<" ";        
        }
        cout<<k<<endl;*/
        int flag=1,tmp=k;
        for(int i=0;i<tmp/2;i++){
            if(arr[i]!=arr[k-1]){
                flag=0;
                break;
            }else{
                k--;
            }
        }
        if(flag==1){
            cout<<"yes !\n";
        }else{
            cout<<"no...\n";
        }
    }
    return 0;
}