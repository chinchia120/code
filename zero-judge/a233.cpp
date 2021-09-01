#include <iostream>
#include <algorithm>
using namespace std;
int main(void){
    int n;
    scanf("%d",&n);
    int num[n];
    for(int i=0;i<n;i++){
        int _num;
        scanf("%d",&_num);
        num[i]=_num;
    }
    sort(num,num+n);
    for(int i=0;i<n;i++){
        printf("%d ",num[i]);
    }
    printf("\n");
}