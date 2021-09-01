#include <bits\stdc++.h>
using namespace std;
bool poker(int *card){
	int a=*card,b=0,c=0,d=0,flag=0;
	for(int i=0;i<5;i++){
		if(*(card+i)==a){
			c++;
		}else{
			if(flag==0){
				b=*(card+i);
				flag=1;
				d++;
				continue;
			}
			if(flag==1){
				if(b==*(card+i)){
					d++;
				}else{
					return false;
				}
			}
		}
	}
	if((c==2 && d==3) || (c==3 && d==2)){
		return true;
	}else{
		return false;
	}
}
int main() {
	int n;
	while(cin>>n){
		for(int i=0;i<n;i++){
			int arr[5]={};
			for(int j=0;j<5;j++){
				cin>>arr[j];
				//cout<<arr[j];
			}
			if(poker(arr)){
				cout<<"Yes\n";
			}else{
				cout<<"No\n"; 
			}
		}
	}
    return 0;
}
