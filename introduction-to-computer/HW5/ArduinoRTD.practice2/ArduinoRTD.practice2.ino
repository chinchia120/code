void setup(){
  pinMode(8,INPUT);
  pinMode(9,INPUT);
  pinMode(10,INPUT);
  pinMode(2,OUTPUT);
}
void loop(){
  int f[3]={262,294,330};
  if(digitalRead(8)==HIGH){
    tone(2,f[0],500);
    delay(1000);
  }else if(digitalRead(9)==HIGH){
    tone(2,f[1],500);
    delay(1000);
  }else if(digitalRead(10)==HIGH){
    tone(2,f[2],500);
    delay(1000);  
  }  
}