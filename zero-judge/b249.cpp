#include <iostream>
using namespace std;

int main(){
    int day;
    while(cin>>day){
        int num[day];
        int sum=0;
        for(int i=0;i<day;i++){
            cin>>num[i];
            for(int j=1;j<day+1;j++){
                sum=sum+num[i]*j;
            }
            cout<<sum<<endl;
        }
    }
    return 0;
}