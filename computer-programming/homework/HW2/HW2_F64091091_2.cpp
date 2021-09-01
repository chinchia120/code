#include <bits\stdc++.h>
#include <windows.h>
void SetColor(int color = 7){
  	HANDLE hConsole;
  	hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
  	SetConsoleTextAttribute(hConsole,color);
}
using namespace std;

int main() {
	string guess,ans="8462";
	int cnt=0;
	while(true){
		cout<<"please enter four different digits number :";
		cin>>guess;
		int A=0,B=0,er=0;
		cnt++;
		if(guess=="0000"){
			return 0;
		}else if(guess=="8462"){
			cout<<"you win in "<<cnt<<" time(s)."<<endl;
			return 0;
		}else if(guess.length()!=4){
			cout<<"you enter "<<guess.length()<<" digit(s) number, please enter ";
			SetColor(14);
			cout<<"\"four\"";
			SetColor();
			cout<<" different digits number."<<endl;
			cout<<"now you guess "<<cnt<<" time(s)"<<endl;
			cout<<endl;
			continue;
		}else{
			for(int i=0;i<4;i++){
				for(int j=i;j<4;j++){
					if(guess[i]==guess[j+1]){
						er=1;
					}
				}
			}
			if(er==1){
				cout<<"you enter repeated number, please enter four ";
				SetColor(14);
				cout<<"\"different\"";
				SetColor();
				cout<<" digits number."<<endl;
				cout<<"now you guess "<<cnt<<" time(s)"<<endl;
				cout<<endl;
				continue;
			}else{
				for(int i=0;i<4;i++){
					for(int j=0;j<4;j++){
						if(i==j && guess[i]==ans[j]){
							A=A+1;
						}else if(i!=j && guess[i]==ans[j]){
							B=B+1;
						}
					}
				}
				cout<<"now you guess "<<cnt<<" time(s) :";
				cout<<A<<"A"<<B<<"B"<<endl;
				cout<<endl;
			}
		}
	}
}
