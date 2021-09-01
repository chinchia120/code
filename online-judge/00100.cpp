#include <iostream>
using namespace std;
int main(void){
    int a,b;
    while(scanf("%d %d",&a,&b)!=EOF){
        printf("%d %d ",a,b);
        if(a>b){
            int tmp=a;
            a=b;
            b=tmp;
        }
        int max=0;
        for(int i=a;i<b+1;i++){
            int cnt=1,num=i;
            while(num!=1){
                if(num%2==0){
                    num=num/2;
                }else{
                    num=num*3+1;
                }
                cnt++;
            }
            if(cnt>max){
                max=cnt;
            }
        }
        printf("%d\n",max);
    }
    return 0;
}