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
        string name[20];
        char gender[20];
        int score[20];
        
    public:
        int cnt=0;
        void Set_Data(string _name,char _gender,int _score){
            name[cnt]=_name;
            gender[cnt]=_gender;
            score[cnt]=_score;
            cnt++;
        }
        void Judge_Data(int n){
            string tmp1;
            char tmp2;
            int tmp3;
            for(int i=0;i<n;i++){
                for(int j=i+1;j<n;j++){
                    if(name[i].length()>name[j].length()){
                        tmp1=name[i];
                        name[i]=name[j];
                        name[j]=tmp1;
                        
                        tmp2=gender[i];
                        gender[i]=gender[j];
                        gender[j]=tmp2;

                        tmp3=score[i];
                        score[i]=score[j];
                        score[j]=tmp3;
                    }else if(name[i].length()==name[j].length()){
                        if(score[i]>score[j]){
                            tmp1=name[i];
                            name[i]=name[j];
                            name[j]=tmp1;

                            tmp2=gender[i];
                            gender[i]=gender[j];
                            gender[j]=tmp2;

                            tmp3=score[i];
                            score[i]=score[j];
                            score[j]=tmp3;
                        }
                    }
                }
            }
        }
        void Output_Data(int n){
            for(int i=0;i<n;i++){
                cout<<name[i]<<" "<<gender[i]<<" "<<score[i]<<endl;
            }
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
        }
        stu.Judge_Data(n);
        stu.Output_Data(n); 
    }
    return 0;
}