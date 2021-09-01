#include <iostream>
#include <string>
using namespace std;
int main(void){
    int n;
    while(cin>>n){
        for(int i=0;i<n;i++){
            string net[10];
            int relate[10],max=0;
            for(int j=0;j<10;j++){
                string _net;
                int _relate;
                cin>>_net>>_relate;
                net[j]=_net;
                relate[j]=_relate;
                if(_relate>max){
                    max=_relate;
                }
            }
            printf("Case #%d:\n",i+1);
            for(int j=0;j<10;j++){
                if(relate[j]==max){
                    cout<<net[j]<<endl;
                }
            }
        }
    }
    return 0;
}
