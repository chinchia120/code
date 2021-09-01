#include <iostream>
#include <iomanip>
using namespace std;

int main() {
	double a,b;
	while(cin>>a){
		b=(a-32)/1.8;
		cout<<fixed<<setprecision(3)<<b<<endl;
	}
    return 0;
    }

