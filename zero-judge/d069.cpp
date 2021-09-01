#include <iostream>
using namespace std;

int main() {
	int n;
	while(cin>>n){
		for(int i=0;i<n;i++){
			int y;
			cin>>y;
			if((y%400==0)or((y%4==0)and(y%100!=0))){
				cout<<"a leap year"<<endl;
			}else{
				cout<<"a normal year"<<endl;
			}
		}
	}
    return 0;
    }

