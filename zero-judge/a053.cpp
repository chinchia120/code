#include <iostream>
using namespace std;

int main() {
	int a,b,c,d;
	cin>>a;
	switch(a){
	case 0 ... 10:
		d=a*6;
		cout<<d<<endl;
		break;
	case 11 ... 19:
		b=a%10;
		d=60+b*2;
		cout<<d<<endl;
		break;
	case 20 ... 39:
		c=a%20;
		d=80+c;
		cout<<d<<endl;
		break;
	case 40 ... 100:
		cout<<100<<endl;
		break;
	}


    return 0;
    }

