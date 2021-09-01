void setup(){
  Serial.begin(9600);
  pinMode(11,OUTPUT);
}
void loop(){
   int value=analogRead(A0);
   Serial.println(value);
   if(value<50 && value>0){
    digitalWrite(11,HIGH);
    delay(1000);
   }else{
    digitalWrite(11,LOW);
    delay(1000);
   }
}