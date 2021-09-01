#include <iostream>
using namespace std;

int main() {
	//question 1
	float x,y;
	cin>>x>>y;
	if(x==0 && y==0){
		cout<<"原點"<<endl;
	}else if(x==0){
		cout<<"在y軸"<<endl;
	}else if(y==0){
		cout<<"在x軸"<<endl; 
	}else if(x>0 && y>0){
		cout<<"第一象限"<<endl;
	}else if(x<0 && y>0){
		cout<<"第二象限"<<endl;
	}else if(x<0 && y<0){
		cout<<"第三象限"<<endl;
	}else if(x>0 && y<0){
		cout<<"第四象限"<<endl;
	}
    return 0;
}

