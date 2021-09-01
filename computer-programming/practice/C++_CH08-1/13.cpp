#include <iostream>
using namespace std;

int main() {
	int m,n,a;
	cin>>m>>n;
	for(int i=0;i<m;i++){
		int max=0,min=101,c=0,d=0;
		for(int j=0;j<n;j++){
			cin>>a;
			if(a>max){
				max=a;
				c=j+1;
			}
			if(a<min){
				min=a;
				d=j+1;
			}
		}
		cout<<d<<" "<<c<<endl;
	}
    return 0;
}
