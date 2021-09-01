#include <iostream>
using namespace std;

int main() {
	int m,n;
	cin>>m>>n;
	int arr[m][n];
	int sum[n]={0};
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			cin>>arr[i][j];
			sum[j]=sum[j]+arr[i][j];
		}
	}
	for(int i=0;i<n;i++){
		cout<<sum[i]<<endl;
	}
    return 0;
}
