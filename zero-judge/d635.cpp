#include <iostream>
#include <vector>
using namespace std;
int main(void){
    int n;
    while(scanf("%d",&n)!=EOF){
        if(n<0){
            cout<<"-1"<<endl;
            return 0;
        }
        vector<int> arr;
        while(n>=8){
            arr.push_back(n%8);
            n/=8;
        }
        cout<<n;
        for(int i=arr.size()-1;i>=0;i--){
            cout<<arr[i];
        }
        cout<<endl;
    }
}