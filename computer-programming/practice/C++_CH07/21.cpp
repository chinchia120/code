#include <bits\stdc++.h>
using namespace std;
void rec(int,int,int &,int &);
int main() {
	int n;
	cin>>n;
	for(int i=0;i<n;i++){
		int l,w,area,peri;
		cin>>l>>w;
		rec(l,w,area,peri);
		cout<<area<<" "<<peri<<endl;
	}
    return 0;
}
void rec (int l,int w,int &area,int &peri){
	area=l*w;
	peri=2*(l+w);
	return;
}

