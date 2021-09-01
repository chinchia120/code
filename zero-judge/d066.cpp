#include <iostream>
using namespace std;

int main() {
	int a,b,c;
	while(cin>>a>>b){
		c=a*60+b;
		switch(c){
			case 450 ... 1019:
				cout<<"At School"<<endl;
				break;
			default:
				cout<<"Off School"<<endl;
			 	break;
		}
	}
    return 0;
    }

