#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;
int main(void){
    int n;
    while(scanf("%d",&n)!=EOF){
        if(n==0){
            return 0;
        }
        vector<int> vec;
        for(int i=0;i<n;i++){
            int a;
            cin>>a;
            vec.push_back(a);
        }
        sort(vec.begin(),vec.end());
        for(int i=0;i<n;i++){
            if(i!=n-1){
                cout<<vec[i]<<" ";
            }else{
                cout<<vec[i]<<endl;
            }
            
        }
    }
    return 0;
}