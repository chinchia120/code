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
    private:
        string name;
        char gender;
        int score;

    public:
        int Lname=0,Sname=1000;
        int Lscore,Sscore;
        void Set_Data(string _name,char _gender,int _score){
            name=_name;
            gender=_gender;
            score=_score;
        }
        void Judge_Data(string name,int score){
            if(name.length()>Lname){
                Lname=name.length();
                Get_Lscore(score);
            }else if(name.length()<Sname){
                Sname=name.length();
                Get_Sscore(score);
            }
        }
        void Get_Lscore(int scoee){
            Lscore=score;
        }
        void Get_Sscore(int score){
            Sscore=score;
        }
        void Output_Data(){
            cout<<Lscore<<endl;
            cout<<Sscore<<endl;
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
            stu.Judge_Data(name,score);
        }
        stu.Output_Data();
    }
    return 0;
}