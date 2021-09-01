#include <bits\stdc++.h>
using namespace std;
struct data{
	string name[20];
	string sex[20];
	int score[20];
}student;
int main(void) {
	int n;
	while(cin>>n){
		string Fname[n],Fsex[n],Mname[n],Msex[n];
		int Fscore[n],Mscore[n];
		int a=0,b=0;
		for(int i=0;i<n;i++){
			cin>>student.name[i]>>student.sex[i]>>student.score[i];
			if(student.sex[i]=="F"){
				Fname[a]=student.name[i];
				Fsex[a]=student.sex[i];
				Fscore[a]=student.score[i];
				a++;
			}else{
				Mname[b]=student.name[i];
				Msex[b]=student.sex[i];
				Mscore[b]=student.score[i];
				b++;
			}
		}
		for(int i=0;i<a;i++){
			for(int j=i+1;j<a;j++){
				if(Fscore[i]<Fscore[j]){
					string tmp1=Fname[i];
					Fname[i]=Fname[j];
					Fname[j]=tmp1;
					
					string tmp2=Fsex[i];
					Fsex[i]=Fsex[j];
					Fsex[j]=tmp2;
					
					int tmp3=Fscore[i];
					Fscore[i]=Fscore[j];
					Fscore[j]=tmp3;
				}
			}
		}
		for(int i=0;i<b;i++){
			for(int j=i+1;j<b;j++){
				if(Mscore[i]<Mscore[j]){
					string tmp4=Mname[i];
					Mname[i]=Mname[j];
					Mname[j]=tmp4;
					
					string tmp5=Msex[i];
					Msex[i]=Msex[j];
					Msex[j]=tmp5;
					
					int tmp6=Mscore[i];
					Mscore[i]=Mscore[j];
					Mscore[j]=tmp6;
				}
			}
		}
		for(int i=0;i<a;i++){
			cout<<Fname[i]<<" "<<Fsex[i]<<" "<<Fscore[i]<<endl;
		}
		for(int i=0;i<b;i++){
			cout<<Mname[i]<<" "<<Msex[i]<<" "<<Mscore[i]<<endl;
		}
	}
    return 0;
}
