#include <iostream>
using namespace std;

int main() {
	int a;
	while(cin>>a){
		switch(a){
			case 0 ... 5:
				cout<<"0"<<endl;
				break;
			case 6 ... 11:
				cout<<"590"<<endl;
				break;
			case 12 ... 17:
				cout<<"790"<<endl;
				break;
			case 18 ... 59:
			    cout<<"890"<<endl;
			    break;
			default:
			    cout<<"399"<<endl;	
				break;				
		}
	}
    return 0;
    }

