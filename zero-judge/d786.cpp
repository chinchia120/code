#include <iostream>
using namespace std;
int main(void){
    int n;
    cin>>n;
    for(int i=0;i<n;i++){
        int m;
        cin>>m;
        double sum=0;
        for(int j=0;j<m;j++){
            int _num;
            cin>>_num;
            sum+=_num;
        }
        printf("%.2lf\n",sum/m);
    }
}
