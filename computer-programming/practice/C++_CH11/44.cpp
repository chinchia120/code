#include <bits\stdc++.h>
using namespace std;
struct data{
	string name[20];
	string sex[20];
	int score[20];
}student;
int main() {
	int n;
	while(cin>>n){
		for(int i=0;i<n;i++){
			cin>>student.name[i]>>student.sex[i]>>student.score[i];
		}
		for(int i=0;i<n;i++){
			for(int j=i+1;j<n;j++){
				if(student.score[i]<student.score[j]){
					int tmp1=student.score[i];
					student.score[i]=student.score[j];
					student.score[j]=tmp1;
					
					string tmp2=student.sex[i];
					student.sex[i]=student.sex[j];
					student.sex[j]=tmp2;
					
					string tmp3=student.name[i];
					student.name[i]=student.name[j];
					student.name[j]=tmp3;
				}
			}
		}
		for(int i=0;i<n;i++){
			cout<<student.name[i]<<" "<<student.sex[i]<<" "<<student.score[i]<<endl;;
		}
	}
    return 0;
}

