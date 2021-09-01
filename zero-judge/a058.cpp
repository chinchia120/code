#include <iostream>
using namespace std;

int main() {
	int n;
	int a=0,b=0,c=0;
	while(cin>>n){
		int array[n];
		for(int i=0;i<n;i++){
			cin>>array[i];
			if(array[i]%3==0){
				a=a+1;
			}else if(array[i]%3==1){
				b=b+1;
			}else if(array[i]%3==2){
				c=c+1;
			}
		}
		cout<<a<<" "<<b<<" "<<c<<endl;
	}
    return 0;
}

