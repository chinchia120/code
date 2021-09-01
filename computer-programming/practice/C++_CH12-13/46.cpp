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
class Student{
	private :
		string name;
		char gender;
		int score;	
		
	public :
		int Fmin=100,Mmin=100;
		string Fname,Mname;
		
		void Set_Data(string _name,char _gender,int _score){
			name=_name;
			gender=_gender;
			score=_score;
		}
		void Judge_Data(string name,char gender,int score){
			if(gender=='F'){
				if(score<Fmin){
					Fmin=score;
					Get_Fname(name);
				}
			}else{
				if(score<Mmin){
					Mmin=score;
					Get_Mname(name);
				}
			}
		}
		void Get_Fname(string name){
			Fname=name;
		}
		void Get_Mname(string name){
			Mname=name;
		}
		void Get_name(){
			cout<<Fname<<endl;
			cout<<Mname<<endl;
		}
};
int main(void){
	Student stu;
	int n;
    while(scanf("%d",&n)!=EOF){
        for(int i=0;i<n;i++){
		    string name;
		    char gender;
		    int score;
		    cin>>name>>gender>>score;
		    stu.Set_Data(name,gender,score);
		    stu.Judge_Data(name,gender,score);
	    }
        stu.Get_name();
    }
    return 0;
}