#include <iostream>
using namespace std;

int main() {
	//question 2
	int time;
	cin>>time;
	if(time<=50){
		cout<<time*100<<endl;
	}else if(time<70){
		cout<<5000+(time-50)*125<<endl;
	}else{
		cout<<7500+(time-70)*150<<endl;
	}
    return 0;
}

