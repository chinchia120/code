#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;
int count(int a,int b){
    if(a>b){
        return a-b;
    }else{
        return b-a;
    }
}
int main(void){
    int n;
    cin>>n;
    for(int i=0;i<n;i++){
        int m;
        cin>>m;
        vector<int> num;
        for(int j=0;j<m;j++){
            int _num;
            cin>>_num;
            num.push_back(_num);
        }
        sort(num.begin(),num.end());
        int mid=num[m/2],sum=0;
        for(int j=0;j<num.size();j++){
            int _sum=count(mid,num[j]);
            sum+=_sum;
        }
        cout<<sum<<endl;
    }
    return 0;
}