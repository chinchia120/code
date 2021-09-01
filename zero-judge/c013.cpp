#include <iostream>
using namespace std;
int main(void){
    int n;
    while(cin>>n){
        for(int i=0;i<n;i++){
            int amp=0,fre=0;
            cin>>amp>>fre;
            for(int j=0;j<fre;j++){
                for(int k=1;k<=amp;k++){
                	for(int l=1;l<=k;l++){
                		cout<<k;
					}
					cout<<endl;
				}
				for(int k=amp-1;k>=1;k--){
                	for(int l=0;l<=k-1;l++){
                		cout<<k;
					}
					cout<<endl;
				}
				cout<<endl;
            }
        }
    }
    return 0;
}