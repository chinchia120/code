#include <iostream>
using namespace std;
int main(void){
    int day;
    while(cin>>day){
        int num[day],sum=0;
        for(int i=0;i<day;i++){
            cin>>num[i];
        }
        for(int i=0;i<day;i++){
            sum=sum+num[i]*(i+1);
        }
        cout<<sum<<endl;
    }
    return 0;
}