#include <iostream>
using namespace std;
int main(void){
    int n;
    while(scanf("%d",&n)!=EOF){
        int arr[n],rev[n-1],cnt=n-1,tmp1=0,tmp2=n,tmp3=0,tmp4=n-1;
        for(int i=0;i<n;i++){
            cin>>arr[i];
            //cout<<arr[i]<<" ";
        }
        //cout<<endl;
        for(int i=0;i<n-1;i++,cnt--){
            rev[i]=arr[cnt];
            //cout<<rev[i]<<" ";
        }
        //cout<<endl;
        for(int i=0;i<n;i++){
            if(i%2==0){
                for(int j=tmp1;j<tmp2;j++){
                    cout<<arr[j]<<" ";
                }
                cout<<endl;
                tmp1++;
                tmp2--;
            }else{
                for(int j=tmp3;j<tmp4;j++){
                	cout<<rev[j]<<" ";
				}
				cout<<endl;
				tmp3++;
                tmp4--;
            }
        }
    }
    return 0;
}