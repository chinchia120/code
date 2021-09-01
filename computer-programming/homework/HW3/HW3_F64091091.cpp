#include <iostream>
#include <ctype.h>
#include <string.h>
#define lenp p.length()
#define lenq q.length()
using namespace std;

int main() {
	while(1){
		string p,q;
		cout<<"Please Enter String p=";
		cin>>p;
		cout<<"Please Enter String q=";
		cin>>q;
		if(p=="0" && q=="0"){
			cout<<"End of Program."<<endl;
			return 0;
		}
		int flag=0;
		for(int i=0;i<lenp;i++){
			if(p[i]<48){
				flag=1;
				break;
			}else if(p[i]>57 && p[i]<65){
				flag=1;
				break;
			}else if(p[i]>90 && p[i]<97){
				flag=1;
				break;
			}else if(p[i]>122){
				flag=1;
				break;
			}
		}
		for(int i=0;i<lenq;i++){
			if(q[i]<48){
				flag=1;
				break;
			}else if(q[i]>57 && q[i]<65){
				flag=1;
				break;
			}else if(q[i]>90 && q[i]<97){
				flag=1;
				break;
			}else if(q[i]>122){
				flag=1;
				break;
			}
		}
		if(flag==1){
			cout<<"Please Enter English Alphabet or Number.\n\n";
			continue;
		}
		//question A
		cout<<"A. ";
		cout<<p<<endl;
		//question B
		cout<<"B. ";
		flag=0;
		for(int i=0;i<lenp;i++){
			if(isupper(p.at(i))){
				cout<<p[i];
				flag=1;
			}
		}
		if(flag==0){
			cout<<"none";
		}
		cout<<"\n";
		//question C
		cout<<"C. ";
		cout<<p+q<<endl;
		//question D 
		cout<<"D. ";
		for(int i=lenp+lenq-1;i>=0;i--){
			cout<<(p+q).at(i);
		}
		cout<<"\n";
		//question E
		cout<<"E. ";
		int sum=0;
		for(int i=0;i<lenp+lenq;i++){
			if(islower((p+q).at(i))){ 
			sum+=1;
			} 
		}
		printf("%d\n",sum);
		//question F
		cout<<"F. ";
		int total=0;
		for(int i=0;i<lenp+lenq;i++){
			if(isdigit((p+q).at(i))){
				total=total+(p+q).at(i)-'0';
			}
		}
		printf("%d\n",total);
		//question G
		cout<<"G. ";
		flag=0;
		for(int i=0;i<lenp;i++){
			for(int j=0;j<lenq;j++){
				if(p.at(i)==q.at(j)){
					cout<<p.at(i);
					flag=1;
				}
			}
		}
		if(flag==0){
			cout<<"none";
		}
		cout<<"\n";
		//question H
		cout<<"H. ";
		int arr1[lenp+lenq]={},arr2[lenq]={};
		for(int i=0;i<lenp;i++){
			arr1[i]=p.at(i);
		}
		for(int i=0;i<lenq;i++){
			arr2[i]=q.at(i);
		}
		int k=0;
		for(int i=0;i<lenq;i++){
			flag=1;
			for(int j=0;j<lenp+lenq;j++){
				if(arr2[i]==arr1[j]){
					flag=0;
					break;
				}
			}
			if(flag==1){
				arr1[lenp+k]=arr2[i];
				k++;
			}
		}
		for(int i=0;i<lenp+lenq;i++){
			cout<<(char)arr1[i];
		}
		cout<<"\n";
		//question I
		cout<<"I. ";
		for(int i=0;i<lenq;i++){
			if(isupper(q.at(i))){
				cout<<char(q.at(i)+32);
			}else if(islower(q.at(i))){
				cout<<char(q.at(i)-32);
			}else{
				cout<<q.at(i);
			}
		}
		cout<<"\n";
		//question J
		cout<<"J. ";
		for(int i=0;i<lenp;i++){
			for(int j=i+1;j<lenp;j++){
				if(p[i]>p[j]){
					char tmp=p[i];
					p[i]=p[j];
					p[j]=tmp;
				}
			}
		}
		for(int i=0;i<lenp;i++){
			cout<<p[i];
		}
		//
		cout<<"\n"<<"\n";
		cout<<"If You Want to End This Program , Please Enter p=0 and q=0.";
		cout<<"\n";
	}
}
