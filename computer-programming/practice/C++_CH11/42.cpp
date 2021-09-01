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
		int Fcnt=0,Mcnt=0;
		for(int i=0;i<n;i++){
			cin>>student.name>>student.sex>>student.score;
			if(student.sex=="F"){
				if(student.score>59){
					Fcnt++;
				}
			}else{
				if(student.score>59){
					Mcnt++;
				}
			}	
		}
		cout<<Fcnt<<endl;
		cout<<Mcnt<<endl;
	}
    return 0;
}

