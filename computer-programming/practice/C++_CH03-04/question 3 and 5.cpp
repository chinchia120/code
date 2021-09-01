#include <iostream>
using namespace std;

int main() {
	//question 3 and 5
	int h,i,j;
	cout<<"input a positive integer A=";
	cin>>h;
	cout<<"input a positive integer B=";
	cin>>i;
	switch(h%i){
		case 0:
			cout<<h/i<<endl;
			break;
		default:
			cout<<h/i+1<<endl;
			break;		
	}
    return 0;
}

