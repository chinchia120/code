#include <iostream>
using namespace std;

int main() {
	int num,range;
	while(cin>>num>>range){
		int sum,time;
		sum=num;
		time=1;
		while(sum<=range){
			sum=sum+num+time;
			time=time+1;
		}
		cout<<time<<endl;
	}
    return 0;
}

