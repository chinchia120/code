#include <iostream>
using namespace std;
int main(void){
    int a,b;
    while(scanf("%d %d",&a,&b)!=EOF){
    	int cnt=0;
    	if(a>b){
    		cout<<"0"<<endl;
    		continue;
		}
        for(int i=a;i<=b;i++){
            int flag=1; 
            for(int j=2;j<i;j++){
                if(i%j==0){
                    flag=0;
                    break;
                }
            }
            if(flag==1){
                cnt++;
            }
        }
        if(a==1){
        	printf("%d\n",cnt-1);
		}else{
			printf("%d\n",cnt);
		}
    }
    return 0;
}