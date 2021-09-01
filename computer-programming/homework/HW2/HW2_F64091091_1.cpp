#include <bits\stdc++.h>
#include <windows.h>
void SetColor(int color = 7){
  	HANDLE hConsole;
  	hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
  	SetConsoleTextAttribute(hConsole,color);
}
using namespace std;

int main() {
	int n,cnt=0;
	while(true){
		cout<<"please enter a number :";
		cin>>n;
		cnt++;
		if(n==66){
			cout<<"you win in "<<cnt<<" times."<<endl;
			return 0;
		}else if(n>66){
			cout<<n<<" is ";
			SetColor(14);
			cout<<"\"bigger\"";
			SetColor();
			cout<<" than the answer."<<endl;
			cout<<"now you guess "<<cnt<<" time(s)"<<endl;
			cout<<endl;
		}else if(n<66){
			cout<<n<<" is ";
			SetColor(14);
			cout<<"\"smaller\"";
			SetColor();
			cout<<" than the answer."<<endl;
			cout<<"now you guess "<<cnt<<" time(s)"<<endl;
			cout<<endl;
		}
	}
}
