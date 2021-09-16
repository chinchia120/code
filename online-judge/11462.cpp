#include <iostream>
#include <vector>
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
        sort(vec,n);
        for(int i=0;i<n;i++){
            cout<<vec[i]<<" ";
        }
        cout<<endl;
    }
}