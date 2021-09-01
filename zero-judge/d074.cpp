#include <iostream>
using namespace std;

int main() {
	int n;
	while(cin>>n){
		int c=0,max=0;
		for(int i=0;i<n;i++){
			cin>>c;
			if(c>max){
				max=c;
			}
		}
		cout<<max<<endl;
	}
    return 0;
    }

