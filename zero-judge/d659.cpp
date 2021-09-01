#include <iostream>
using namespace std;
int main(void){
    int n;
    while(scanf("%d",&n)!=EOF){
        for(int i=0;i<n;i++){
            int num[3]={};
            for(int j=0;j<3;j++){
                int _num;
                cin>>_num;
                num[j]=_num;
            }
            for(int j=0;j<3;j++){
                for(int k=j+1;k<3;k++){
                    if(num[j]>num[k]){
                        int tmp=num[j];
                        num[j]=num[k];
                        num[k]=tmp;
                    }
                }
            }
            printf("Case %d: %d\n",i+1,num[1]);
        }
    }
}