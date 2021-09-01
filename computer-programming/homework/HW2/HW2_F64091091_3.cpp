#include <iostream>
#include <time.h>
#include <algorithm>
#include <windows.h>
#include <cstdlib>
#include <string> 
void SetColor(int color = 7){
  	HANDLE hConsole;
  	hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
  	SetConsoleTextAttribute(hConsole,color);
}
using namespace std;
#define rs random_shuffle 
int main() {
	//make a number
	srand(time(0));
	int ans=0,arr[4],cnt=0,tmp;
	while(true){
		int er=0;
		ans=rand()%10000;
	 	tmp=ans;
		for(int i=0;i<4;i++,ans/=10){
			arr[i]=ans%10;
		}
		for(int i=0;i<4;i++){
			for(int j=i;j<4;j++){
				if(arr[i]==arr[j+1]){
					er=1;
				}
			}
		}
		if(er==1){
			continue;
		}else{
			//cout<<tmp<<endl;
			break;
		}
	}
	//guess
	while(true){
		string gu;
		int ans=tmp,ran[4];
		cout<<"please enter four different digits number :";
		cin>>gu;
		int guess=atoi(gu.c_str());
		cnt++;
		if(guess==ans){
			break;
		}else if(guess==0000){
			return 0;
		}else if(gu.length()!=4){
			cout<<"you enter "<<gu.length()<<" digit(s) number, please enter ";
			SetColor(14);
			cout<<"\"four\"";
			SetColor();
			cout<<" different digits number."<<endl;
			cout<<"now you guess "<<cnt<<" time(s)"<<endl;
			cout<<endl;
		}else{
			for(int i=0;i<4;i++,guess=guess/10,ans=ans/10){
				if(guess%10>ans%10){
					ran[i]=guess%10-ans%10;
				}else{
					ran[i]=-(guess%10-ans%10);
				}
			}
			rs(ran,ran+4);
			cout<<"the difference of each digit after arbitrary arrangement are :";
			for(int i=0;i<4;i++){
				cout<<ran[i]<<" \n"[i==3];
			}
			cout<<"now you guess "<<cnt<<" time(s)"<<endl;
			cout<<endl;
		}
	}
	cout<<"you win in "<<cnt<<" time(s)."<<endl;
	return 0;
}
