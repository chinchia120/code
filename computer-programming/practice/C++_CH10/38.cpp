#include <iostream>
using namespace std;
int main(void){
	int n;
	while(cin>>n){
		for(int i=0;i<n;i++){
			string str;
			cin>>str;
			for(int j=0;j<str.length();j++){
				for(int k=j+1;k<str.length();k++){
					if(str[j]>str[k]){
						char tmp=str[j];
						str[j]=str[k];
						str[k]=tmp;
					}
				}
			}
			for(int j=0;j<str.length();j++){
				cout<<str[j];
			}
			cout<<"\n";
		}
	}
	return 0;
}

