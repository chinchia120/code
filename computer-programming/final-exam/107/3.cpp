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
            int a;
            cin>>a;
            for(int j=2;j<a+1;j++){
                int q=a-j,flag1=1,flag2=1;
                for(int k=2;k<j;k++){
                	if(j%k==0){
                		flag1==0;
                		break;
					}
				}
				for(int k=2;k<q;k++){
					if(q%k==0){
						flag2=0;
						break;
					}
				}
				if(flag1==1 && flag2==1){
					cout<<j<<"+"<<q<<endl;
					break;
				}
            }
        }
    }
    return 0;
}