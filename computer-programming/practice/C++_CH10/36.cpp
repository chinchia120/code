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
			int cnt=0;
			for(int j=0;j<str.length();j++){
				if(isalpha(str[j])){
					cnt++;	
				}
			}
			cout<<cnt<<endl;
		}
	}
	return 0;
}

