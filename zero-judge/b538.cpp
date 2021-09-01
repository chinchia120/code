#include <iostream>
using namespace std;
void rf(int a,int b){
    if(b<0){
        a=-a;
        b=-b;
    }
    int tmp=b;
    for(int i=2;i<=tmp;i++){
        while(a%i==0 && b%i==0){
            a=a/i;
            b=b/i;
        }
    }
    if(a==0 || b==1){
        cout<<a<<endl;
    }else{
        cout<<a<<"/"<<b<<endl;
    }
    return;
}
int main(void){
    int a,b,c,d;
    char ch;
    while(cin>>a>>b>>c>>d>>ch){
        int nume,deno;
        if(ch=='+'){
            nume=a*d+b*c;
            deno=b*d;
            rf(nume,deno);
        }else if(ch=='-'){
            nume=a*d-b*c;
            deno=b*d;
            rf(nume,deno);
        }else if(ch=='*'){
            nume=a*c;
            deno=b*d;
            rf(nume,deno);
        }else{
            nume=a*d;
            deno=b*c;
            rf(nume,deno);
        }
    }
    return 0;
}