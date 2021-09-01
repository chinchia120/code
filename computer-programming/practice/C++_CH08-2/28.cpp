#include <bits\stdc++.h>
using namespace std;
int mid3(int arr[],int m){
	for(int i=0;i<m;i++){
		for(int j=i+1;j<m;j++){
			if(arr[i]>arr[j]){
				int tmp=arr[i];
				arr[i]=arr[j];
				arr[j]=tmp;
			}
		}
	}
	return arr[m/2];
}
int main() {
	int n;
	while(cin>>n){
		for(int i=0;i<n;i++){
			int m;
			cin>>m;
			int arr[m]={};
			for(int j=0;j<m;j++){
				cin>>arr[j];
			}
			cout<<mid3(arr,m)<<endl;
		}
	}
    return 0;
}
