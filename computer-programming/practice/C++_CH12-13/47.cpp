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
        int Fcnt=0,Mcnt=0;
        void Set_Data(string _name,char _gender,int _score){
            name=_name;
            gender=_gender;
            score=_score;
        }
        void Judge_Data(char gender,int score){
            if(gender=='F'){
                if(score<60){
                    Fnum();
                }
            }else{
                if(score<60){
                    Mnum();
                }
            }
        }
        void Fnum(){
            Fcnt++;
        }
        void Mnum(){
            Mcnt++;
        }      
        void cnt(){
            cout<<Fcnt<<endl;
            cout<<Mcnt<<endl;
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
            stu.Judge_Data(gender,score);
        }
        stu.cnt();
    }
    return 0;
}