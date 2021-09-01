void setup(){
  Serial.begin(9600);
  pinMode(13,OUTPUT);//red
  pinMode(12,OUTPUT);//yellow
  pinMode(11,OUTPUT);//blue
  pinMode(10,OUTPUT);//voice
  pinMode(A0,INPUT);//water
}

void loop() {
  int value;
  value=analogRead(A0);
  Serial.print(value);
  Serial.print("  ");
  if(value>=520){
    Serial.println("HIGH LEVEL");
    digitalWrite(13,HIGH);
    digitalWrite(12,LOW);
    digitalWrite(11,LOW);
    tone(10,1500,200);
  }else if(value>=120 && value<520){
    Serial.println("MIDDLE LEVEL");
    digitalWrite(13,LOW);
    digitalWrite(12,HIGH);
    digitalWrite(11,LOW);
  }else if(value<120){
    Serial.println("LOW LEEVL");
    digitalWrite(13,LOW);
    digitalWrite(12,LOW);
    digitalWrite(11,HIGH);
  }  
  delay(500); 
}