#include <iostream>
using namespace std;
int main(void){
    int n;
    while(cin>>n){
        for(int i=0;i<n;i++){
            float G,A,H,W,BMR;
            cin>>G>>A>>H>>W;
            if(G==1){
                BMR=(13.7*W)+(5*H)-(6.8*A)+66;
                printf("%.2f\n",BMR);
            }else if(G==0){
                BMR=(9.6*W)+(1.8*H)-(4.7*A)+655;
                printf("%.2f\n",BMR);
            }
        }
    }
    return 0;
}