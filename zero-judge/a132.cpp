#include <iostream>
#include <vector>
using namespace std;
int main(void){
    int n;
    while(cin>>n){
        if(n==0){
            return 0;
        }
        vector<int> arr;
        while(n!=1){
            int tmp=n%2;
            arr.push_back(tmp);
            n/=2;
        }
        vector<int> num;
        num.push_back(1);
        for(int i=arr.size()-1;i>=0;i--){
            num.push_back(arr[i]);
        }
        int cnt=0;
        for(int i=0;i<num.size();i++){
            if(num[i]==1){
                cnt++;
            }
        }
        cout<<"The parity of ";
        for(int i=0;i<num.size();i++){
            cout<<num[i];
        }
        printf(" is %d (mod 2).\n",cnt);
    }
}