#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;
int main(void){
    int n;
    while(cin>>n){
        string arr[n];
        for(int i=0;i<n;i++){
            cin>>arr[i];
            //cout<<arr[i]<<" ";
        }
        int a[n];
        for(int i=0;i<n;i++){
            a[i]=atoi(arr[i].c_str());
            //cout<<a[i]<<" ";
        }
        string str[n];
        for(int i=0;i<n;i++){
            str[i]=arr[i].at(arr[i].length()-1);
            //cout<<str[i]<<" ";
        }
        for(int i=0;i<n;i++){
            for(int j=i+1;j<n;j++){
                if(str[i]>str[j]){
                    int tmp1=a[i];
                    a[i]=a[j];
                    a[j]=tmp1;
                }else if(str[i]==str[j]){
                    if(a[i]<a[j]){
                        int tmp2=a[i];
                        a[i]=a[j];
                        a[j]=tmp2;
                    }
                }
            }
        }
        for(int i=0;i<n;i++){
            cout<<a[i]<<" ";
        }
        cout<<endl;
    }
    return 0;
}