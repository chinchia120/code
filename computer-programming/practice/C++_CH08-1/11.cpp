#include <iostream>
using namespace std;

int main() {
	int n,max=0,min=101;
	cin>>n;
	for(int i=0;i<n;i++){
		int a;
		cin>>a;
		if(a>max){
			max=a;
		}
		if(a<min){
			min=a;
		}
	}
	cout<<min<<" "<<max<<endl;
    return 0;
}
