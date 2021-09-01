#include <iostream>
#include <math.h>
using namespace std;

int main() {
	//question 2 and 4
	float c,d,e;
	int f,g;
	cout<<"input a positive integer A=";
	cin>>c;
	cout<<"input a positive integer B=";
	cin>>d;
	e=c/d;
	f=c/d;
	g=floor((e-f)*10);
	switch(g){
		case 0 ... 4:
			cout<<f<<endl;
			break;
		default:
			cout<<f+1<<endl;
			break; 	
	}	
    return 0;
}

