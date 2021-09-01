#include <bits\stdc++.h>
using namespace std;

int main() {
	int n;
	while(cin>>n){
		for(int i=0;i<n;i++){
			string str;
			cin>>str;
			for(int j=0;j<str.length();j++){
				if(isupper(str.at(j))){
					cout<<char(str.at(j)+32);
				}else{
					cout<<str.at(j);
				}
			}
			cout<<endl;
		}
	}
    return 0;
}
