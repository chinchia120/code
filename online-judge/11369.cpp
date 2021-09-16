#include <iostream>
#include <vector>
#include <algorithm>
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
        sort(vec.begin(),vec.end(),greater<int>());
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