#include <iostream>
using namespace std;
int main(void){
    int n;
    while(cin>>n){
        for(int i=0;i<n;i++){
            int m;
            cin>>m;
            int num[m]={};
            for(int j=0;j<m;j++){
                int _num;
                cin>>_num;
                num[j]=_num;
            }
            int inc=0,dec=0;
            for(int j=0;j<m-1;j++){
                if(num[j+1]-num[j]>0){
                    inc++;
                }else if(num[j+1]-num[j]<0){
                    dec++;
                }
            }
            printf("Case %d: %d %d\n",i+1,inc,dec);
        }
    }
    return 0;
}