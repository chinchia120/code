#include <iostream>
using namespace std;
int main(void){
    int ans=0,n=0;
    scanf("%d %d",&ans,&n);
    int tmp=ans;
    for(int i=0;i<n;i++){
        int guess=0,gu[4];
        scanf("%d",&guess);
        for(int j=0;j<4;j++,guess/=10){
            gu[j]=guess%10;
            //cout<<gu[i]<<" ";
        }
        int ANS[4],ans=tmp;
        for(int i=0;i<4;i++,ans/=10){
            ANS[i]=ans%10;
            //cout<<ANS[i]<<" ";
        }
        int A=0,B=0;
        for(int j=0;j<4;j++){
            for(int k=0;k<4;k++){
                if(j==k && ANS[j]==gu[k]){
                    A=A+1;
                    ANS[j]=10;
                }
            }
        }
        for(int j=0;j<4;j++){
            for(int k=0;k<4;k++){
                if(j!=k && ANS[j]==gu[k]){
                    B=B+1;
                    ANS[j]=10;
                }
            }
        }
        printf("%dA%dB\n",A,B);
    }
    return 0;
}