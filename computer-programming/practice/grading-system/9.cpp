#include <iostream>
using namespace std;

int main() {
	int n,a;
	while(cin>>n){
		for(int i=0;i<n;i++){
			cin>>a;
			for(int j=0;j<a;j++){
				for(int k=1;k<=a-j;k++){
					cout<<"*"; 
				}
				cout<<endl;
			}	
		}
	}
    return 0;
}

