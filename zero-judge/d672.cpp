#include <iostream>
#include <string>
using namespace std;
int count_digit(int sum){
    int cnt=0;
    while(sum!=0){
        sum/=10;
        cnt++;
    }
    return cnt;
}
int count_sum(int sum){
    int _sum=0;
    int digit=count_digit(sum);
    for(int i=0;i<digit;i++){
        _sum+=sum%10;
        sum/=10;
    }
    return _sum;
}
int main(void){
    string str;
    while(cin>>str){
        if(str=="0"){
            return 0;
        }
        if(str=="9"){
            cout<<str<<" is a multiple of 9 and has 9-degree 1."<<endl;
        }else{
            int sum=0;
            for(int i=0;i<str.length();i++){
                sum+=str[i]-'0';
            }
            int cnt=1;
            while(sum>9){
                sum=count_sum(sum);
                cnt++;
            }
            if(sum==9){
                cout<<str<<" is a multiple of 9 and has 9-degree "<<cnt<<"."<<endl;
            }else{
                cout<<str<<" is not a multiple of 9."<<endl;
            }
        }
    }
}