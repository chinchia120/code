#include <iostream>
using namespace std;
int main(void){
    int n;
    cin>>n;
    for(int i=0;i<n;i++){
        int e,f,c;
        cin>>e>>f>>c;
        int tmp=e+f,sum=0;
        while(tmp>=c){
            sum+=tmp/c;
            tmp=tmp/c+tmp%c;
        }
        cout<<sum<<endl;
    }
}