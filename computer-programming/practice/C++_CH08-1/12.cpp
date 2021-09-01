#include <iostream>
using namespace std;

int main() {
	int m,n,a;
	cin>>m>>n;
	for(int i=0;i<m;i++){
		int max=0,min=101;
		for(int j=0;j<n;j++){
			cin>>a;
			if(a>max){
				max=a;
			}
			if(a<min){
				min=a;
			}
		}
		cout<<min<<" "<<max<<endl;
	}
    return 0;
}

