#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <stack>
#include <algorithm>
#include <vector>
using namespace std;

int main() {
    int n;
    while(scanf("%d",&n)!=EOF){
        for(int i=0;i<n;i++){
            int a,b,c;
            char p,q;
            cin>>a>>p>>b>>q>>c;
            if(p=='+'){
                if(q=='+'){
                    cout<<a+b+c<<endl;
                }else if(q=='-'){
                    cout<<a+b-c<<endl;
                }else if(q=='*'){
                    cout<<a+b*c<<endl;
                }else{
                    cout<<a+b/c<<endl;
                }
            }else if(p=='-'){
                if(q=='+'){
                    cout<<a-b+c<<endl;
                }else if(q=='-'){
                    cout<<a-b-c<<endl;
                }else if(q=='*'){
                    cout<<a-b*c<<endl;
                }else{
                    cout<<a-b/c<<endl;
                }
            }else if(p=='*'){
                if(q=='+'){
                    cout<<a*b+c<<endl;
                }else if(q=='-'){
                    cout<<a*b-c<<endl;
                }else if(q=='*'){
                    cout<<a*b*c<<endl;
                }else{
                    cout<<a*b/c<<endl;
                }
            }else{
                if(q=='+'){
                    cout<<a/b+c<<endl;
                }else if(q=='-'){
                    cout<<a/b-c<<endl;
                }else if(q=='*'){
                    cout<<a/b*c<<endl;
                }else{
                    cout<<a/b/c<<endl;
                }
            }
        }
    }
    return 0;
}
