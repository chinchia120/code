#include <iostream>
#include <vector>
using namespace std;
int main(void){
    int n;
    while(cin>>n){
        int cnt=0;
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
        int flag=1;
        for(int i=0;i<n-1;i++){
            if((vec[i+1]-vec[i])!=1){
                flag=0;
                break;
            }
        }
        if(flag==1){
            printf("%d %d yes\n",vec[0],vec[n-1]);
        }else{
            printf("%d %d no\n",vec[0],vec[n-1]);
        }
    }
    return 0;
}