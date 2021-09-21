#include <iostream>
#include <string>
using namespace std;
int main(void){
    string _n;
    while(getline(cin,_n)){
        int n=stoi(_n);
        for(int i=0;i<n;i++){
            string str1,str2;
            getline(cin,str1);
            getline(cin,str2);
            int len1=str1.length();
            int len2=str2.length();
            int flag=1;
            if(len1==len2){
                for(int j=0;j<len1;j++){
                    if(str1[j]!=str2[j]){
                        flag=0;
                        break;
                    }
                }
                if(flag==1){
                    printf("Case %d: Yes\n",i+1);
                }else{
                    printf("Case %d: Wrong Answer\n",i+1);
                }
            }else{
                int cnt=0;
                for(int j=0;j<len1;j++){
                    if(isspace(str1[j])){
                        continue;
                    }else{
                        if(str1[j]!=str2[cnt]){
                            flag=0;
                            break;
                        }else{
                            cnt++;
                        }
                    }
                }
                if(flag==1){
                    printf("Case %d: Output Format Error\n",i+1);
                }else{
                    printf("Case %d: Wrong Answer\n",i+1);
                }
            }
        }
    }    
    return 0;
}