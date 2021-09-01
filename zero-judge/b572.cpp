#include <iostream>
using namespace std;

int main() {
	int n;
	while(cin>>n){
		int a,b,c,d,e,begin=0,leave=0;
		for(int i=0;i<n;i++){
			cin>>a>>b>>c>>d>>e;
			begin=a*60+b;
			leave=c*60+d;
			if(leave-begin>=e){
				cout<<"Yes"<<endl;
			}else{
				cout<<"No"<<endl;
			}
		}
	}
    return 0;
    }

