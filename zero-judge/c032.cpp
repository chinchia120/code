#include <iostream>
using namespace std;
int main(void){
    int ans;
    cout<<"PERFECTION OUTPUT\n";
    while(cin>>ans){
        if(ans==0){
            cout<<"END OF OUTPUT\n";
            return 0;
        }else{
            long long unsigned int sum=0;
            for(int i=1;i<ans;i++){
                if(ans%i==0){
                    sum=sum+i;
                }
            }
            //cout<<sum<<endl;
            if(sum==ans){
                printf("%5d  PERFECT\n",ans);
            }else if(sum<ans){
                printf("%5d  DEFICIENT\n",ans);
            }else{
                printf("%5d  ABUNDANT\n",ans);
            }
        }
    }
}