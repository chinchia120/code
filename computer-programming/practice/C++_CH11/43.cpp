#include <bits\stdc++.h>
using namespace std;
struct data{
	string name;
	string sex;
	int score;
}student;
int main() {
	int n;
	while(cin>>n){
		int Lname=0,Sname=999;
		string Lsex,Ssex;
		for(int i=0;i<n;i++){
			cin>>student.name>>student.sex>>student.score;
			if((student.name).length()>Lname){
				Lname=(student.name).length();
				Lsex=student.sex;
			}else if((student.name).length()<Sname){
				Sname=(student.name).length();
				Ssex=student.sex;
			}
		}
		cout<<Lsex<<endl;
		cout<<Ssex<<endl;
	}
    return 0;
}

