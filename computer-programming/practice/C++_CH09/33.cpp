#include <bits\stdc++.h>
using namespace std;
void circle(int r,double *area,double *peri){
	*area=r*r*3.1416;
	*peri=2*3.1416*r;
}
int main() {
	int n;
	while(cin>>n){
		for(int i=0;i<n;i++){
			int a;
			double area,peri;
			cin>>a;
			circle(a,&area,&peri);
			printf("%.2f %.2f\n",area,peri);
		}
	}
    return 0;
}

