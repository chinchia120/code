#include <bits\stdc++.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <stack>
#include <algorithm>
#include <vector>
using namespace std;
struct data{
	string name;
	string sex;
	int score;
}student;
int main() {
	int n;
	while(cin>>n){
		int Fmax=0,Mmax=0;
		string Fname,Mname;
		for(int i=0;i<n;i++){
			cin>>student.name>>student.sex>>student.score;
			if(student.sex=="F"){
				if(student.score>Fmax){
					Fmax=student.score;
					Fname=student.name;
				}
			}else{
				if(student.score>Mmax){
					Mmax=student.score;
					Mname=student.name;
				}
			}
		}
		cout<<Fname<<endl;
		cout<<Mname<<endl;
	}
    return 0;
}

