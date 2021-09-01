#include <iostream>
using namespace std;
struct data{
    int x[1000];
    int y[1000];
}position;
int main(void){
    int n;
    while(cin>>n){
        for(int i=0;i<n;i++){
            cin>>position.x[i]>>position.y[i];
        }
        for(int i=0;i<n;i++){
            for(int j=i+1;j<n;j++){
                if(position.x[i]>position.x[j]){
                    int tmp1=position.x[i];
                    position.x[i]=position.x[j];
                    position.x[j]=tmp1;

                    int tmp2=position.y[i];
                    position.y[i]=position.y[j];
                    position.y[j]=tmp2;
                }else if(position.x[i]==position.x[j]){
                    if(position.y[i]>position.y[j]){
                        int tmp3=position.x[i];
                        position.x[i]=position.x[j];
                        position.x[j]=tmp3;

                        int tmp4=position.y[i];
                        position.y[i]=position.y[j];
                        position.y[j]=tmp4;
                    }
                }
            }
        }
        for(int i=0;i<n;i++){
            cout<<position.x[i]<<" "<<position.y[i]<<endl;
        }
    }
    return 0;
}