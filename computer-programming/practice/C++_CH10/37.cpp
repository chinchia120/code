#include <iostream>
#include <string.h>
using namespace std;
int main(void){
	int n;
	while(cin>>n){
		for(int i=0;i<n;i++){
			string str;
			char *ptr;
			ptr=new char[100];
			cin>>str;
			strcpy(ptr,str.c_str());
			for(int j=str.length()-1;j>=0;j--){
				cout<<ptr[j];
			}
			cout<<endl;
		}
	}
	return 0;
}

