#include <iostream>
using namespace std;

int main() {
	//question 4
	string str;
	cin>>str;
	if(str.at(1)=='0' && str.at(2)=='0' && str.at(3)=='0'){
		cout<<str.at(0)<<"�d"<<endl; 
	}else if(str.at(1)=='0' && str.at(2)=='0'){
		cout<<str.at(0)<<"�d�s"<<str.at(3)<<endl; 
	}else if(str.at(2)=='0' && str.at(3)!='0'){
		cout<<str.at(0)<<"�d"<<str.at(1)<<"�ʹs"<<str.at(3)<<endl; 
	}else if(str.at(1)=='0' && str.at(2)!='0' && str.at(3)!='0'){
		cout<<str.at(0)<<"�d�s"<<str.at(2)<<"�Q"<<str.at(3)<<endl;	
	}else if(str.at(2)=='0' && str.at(3)=='0'){
		cout<<str.at(0)<<"�d"<<str.at(1)<<"��"<<endl; 
	}else if(str.at(1)=='0' && str.at(3)=='0'){
		cout<<str.at(0)<<"�d�s"<<str.at(2)<<"�Q"<<endl; 
	}else if(str.at(3)=='0'){
		cout<<str.at(0)<<"�d"<<str.at(1)<<"��"<<str.at(2)<<"�Q"<<endl;
	}else{
		cout<<str.at(0)<<"�d"<<str.at(1)<<"��"<<str.at(2)<<"�Q"<<str.at(3)<<endl; 
	}
    return 0;
}

