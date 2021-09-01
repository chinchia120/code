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
        for(int i=0;i<n;i++){
            for(int j=i+1;j<n;j++){
                if(vec[i]>vec[j]){
                    int tmp=vec[i];
                    vec[i]=vec[j];
                    vec[j]=tmp;
                }
            }
        }
        for(int i=0;i<n;i++){
            cout<<vec[i]<<" ";
        }
        cout<<endl;
    }
}