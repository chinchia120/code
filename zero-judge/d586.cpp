#include <iostream>
#include <ctype.h>
using namespace std;

int main(){
    int n;
    while(cin>>n){
        for(int i=0;i<n;i++){
            int m,sum=0,flag=0;
            cin>>m;
            for(int j=0;j<m;j++){
                string str;
                cin>>str;
                if(str.at(0)>=97 && str.at(0)<=122){
                	flag=1;
					if(str=="u"){
						sum=sum+1;
					}else if(str=="z"){
						sum=sum+2;
					}else if(str=="r"){
						sum=sum+3;
					}else if(str=="m"){
						sum=sum+4;
					}else if(str=="a"){
						sum=sum+5;
					}else if(str=="t"){
						sum=sum+6;
					}else if(str=="i"){
						sum=sum+7;
					}else if(str=="f"){
						sum=sum+8;
					}else if(str=="x"){
						sum=sum+9;
					}else if(str=="o"){
						sum=sum+10;
					}else if(str=="p"){
						sum=sum+11;
					}else if(str=="n"){
						sum=sum+12;
					}else if(str=="h"){
						sum=sum+13;
					}else if(str=="w"){
						sum=sum+14;
					}else if(str=="v"){
						sum=sum+15;
					}else if(str=="b"){
						sum=sum+16;
					}else if(str=="s"){
						sum=sum+17;
					}else if(str=="l"){
						sum=sum+18;
					}else if(str=="e"){
						sum=sum+19;
					}else if(str=="k"){
						sum=sum+20;
					}else if(str=="y"){
						sum=sum+21;
					}else if(str=="c"){
						sum=sum+22;
					}else if(str=="q"){
						sum=sum+23;
					}else if(str=="j"){
						sum=sum+24;
					}else if(str=="g"){
						sum=sum+25;
					}else if(str=="d"){
						sum=sum+26;
					}
                }else{
                    if(str=="1"){
                        cout<<"m";
                    }else if(str=="2"){
                        cout<<"j";
                    }else if(str=="3"){
                        cout<<"q";
                    }else if(str=="4"){
                        cout<<"h";
                    }else if(str=="5"){
                        cout<<"o";
                    }else if(str=="6"){
                        cout<<"f";
                    }else if(str=="7"){
                        cout<<"a";
                    }else if(str=="8"){
                        cout<<"w";
                    }else if(str=="9"){
                        cout<<"c";
                    }else if(str=="10"){
                        cout<<"p";
                    }else if(str=="11"){
                        cout<<"n";
                    }else if(str=="12"){
                        cout<<"s";
                    }else if(str=="13"){
                        cout<<"e";
                    }else if(str=="14"){
                        cout<<"x";
                    }else if(str=="15"){
                        cout<<"d";
                    }else if(str=="16"){
                        cout<<"k";
                    }else if(str=="17"){
                        cout<<"v";
                    }else if(str=="18"){
                        cout<<"g";
                    }else if(str=="19"){
                        cout<<"t";
                    }else if(str=="20"){
                        cout<<"z";
                    }else if(str=="21"){
                        cout<<"b";
                    }else if(str=="22"){
                        cout<<"l";
                    }else if(str=="23"){
                        cout<<"r";
                    }else if(str=="24"){
                        cout<<"y";
                    }else if(str=="25"){
                        cout<<"u";
                    }else if(str=="26"){
                        cout<<"i";
                    }
                }
            }
            if(flag==1){
            	cout<<sum;
			}
			cout<<endl;
        }
    }
    return 0;
}