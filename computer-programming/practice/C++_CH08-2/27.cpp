#include <bits\stdc++.h>
using namespace std;
int mid2(int arr[]){
	for(int j=0;j<5;j++){
		for(int k=j+1;k<5;k++){
			if(arr[j]>arr[k]){
				int tmp=arr[j];
				arr[j]=arr[k];
				arr[k]=tmp;
			}
		}
	}
	return arr[2];
}
int main() {
	int n;
	while(cin>>n){
		int arr[5]={};
		for(int i=0;i<n;i++){
			for(int j=0;j<5;j++){
				cin>>arr[j];
			}
			cout<<mid2(arr)<<endl;
		}
	}
    return 0;
}

