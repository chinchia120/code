#include <iostream>
#include <vector>
using namespace std;
int main(void){
    int n;
    cin>>n;
    for(int i=0;i<n;i++){
        int m;
        cin>>m;
        vector<int> vec;
        for(int j=0;j<m;j++){
            int _num;
            cin>>_num;
            vec.push_back(_num);
        }
        for(int j=0;j<m;j++){
            for(int k=j+1;k<m;k++){
                if(vec[j]<vec[k]){
                    int tmp=vec[j];
                    vec[j]=vec[k];
                    vec[k]=tmp;
                }
            }
        }
        int sum=0;
        for(int j=0;j<m;j++){
            if(j%3==2){
                sum+=vec[j];
            }
        }
        cout<<sum<<endl;
    }
    return 0;
}