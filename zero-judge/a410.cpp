#include <iostream>
#include <iomanip>
using namespace std;

int main() {
	double a,b,c,d,e,f,g,h,i;
	cin>>a>>b>>c>>d>>e>>f;
	g=a*e-b*d;
	h=c*e-b*f;
	i=a*f-c*d;
	if(g==0){
		if(h==0&&i==0){
			cout<<"Too many"<<endl;
		}else{
			cout<<"No answer"<<endl;
		}		
	}else{
		cout<<"x="<<fixed<<setprecision(2)<<h/g<<endl;
		cout<<"y="<<fixed<<setprecision(2)<<i/g<<endl;
	}
    return 0;
    }

