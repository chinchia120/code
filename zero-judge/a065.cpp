#include <iostream>
using namespace std;

int main() {
	string c;
	char ch[7];
	int num[7];
	while(cin>>c){
		for(int i=0;i<7;i++){
			ch[i]=c.at(i);
			num[i]=(int)(ch[i]);
		}
		if(num[1]-num[0]>0 || num[1]-num[0]==0){
			cout<<num[1]-num[0];
		}else if(num[1]-num[0]<0){
			cout<<-(num[1]-num[0]);
		}
		if(num[2]-num[1]>0 || num[2]-num[1]==0){
			cout<<num[2]-num[1];
		}else if(num[2]-num[1]<0){
			cout<<-(num[2]-num[1]);
		}
		if(num[3]-num[2]>0 || num[3]-num[2]==0){
			cout<<num[3]-num[2];
		}else if(num[3]-num[2]<0){
			cout<<-(num[3]-num[2]);
		}
		if(num[4]-num[3]>0 || num[4]-num[3]==0){
			cout<<num[4]-num[3];
		}else if(num[4]-num[3]<0){
			cout<<-(num[4]-num[3]);
		}
		if(num[5]-num[4]>0 || num[5]-num[4]==0){
			cout<<num[5]-num[4];
		}else if(num[5]-num[4]<0){
			cout<<-(num[5]-num[4]);
		}
		if(num[6]-num[5]>0 || num[6]-num[5]==0){
			cout<<num[6]-num[5];
		}else if(num[6]-num[5]<0){
			cout<<-(num[6]-num[5]);
		}
		cout<<endl;
	}
	return 0;
}

