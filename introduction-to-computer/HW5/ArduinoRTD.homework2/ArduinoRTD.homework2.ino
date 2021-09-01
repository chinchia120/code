int f[63]={330,330,330,0,330,330,330,0,330,392,262,294,330,0,0,0,349,349,349,349,349,330,330,330,330,294,294,330,294,0,392,0,330,330,330,0,330,330,330,0,330,392,262,294,330,0,0,0,349,349,349,349,349,330,330,330,392,392,349,294,262,0,0};

void setup(){
  Serial.begin(9600);
  pinMode(13,OUTPUT);
  pinMode(12,INPUT);
  pinMode(11,OUTPUT);
  pinMode(10,INPUT);
  pinMode(9,OUTPUT);
  pinMode(8,INPUT);
  pinMode(7,OUTPUT);
  pinMode(6,INPUT);
  pinMode(5,OUTPUT);
  pinMode(4,INPUT);
  pinMode(3,OUTPUT);
  pinMode(2,OUTPUT);
  pinMode(A5,INPUT);
}

void loop(){
  for(int i=0;i<999999;i++){
    for(int j=0;j<63;j++){
      float duration,distance;
      digitalWrite(2,HIGH);
      delayMicroseconds(10);
      digitalWrite(2,LOW);
      duration=pulseIn(A5,HIGH); 
      distance=(duration/2)*0.034;
      Serial.print(distance);
      Serial.println(" cm");
      delay(500);
      tone(3,f[j],distance*7.5);
      delay(distance*7.5);
      if(digitalRead(12)==HIGH){
        digitalWrite(13,HIGH);
        break;
      }else{
        digitalWrite(13,LOW);
      }
      if(digitalRead(10)==HIGH){
        digitalWrite(11,HIGH);
      }else{
        digitalWrite(11,LOW);
      }
      if(digitalRead(8)==HIGH){
        digitalWrite(9,HIGH);
      }else{
        digitalWrite(9,LOW);
      }
      if(digitalRead(6)==HIGH){
        digitalWrite(7,HIGH);
      }else{
        digitalWrite(7,LOW);
      }
      if(digitalRead(4)==HIGH){
      digitalWrite(5,HIGH);
      }else{
        digitalWrite(5,LOW);
      }
    }
    if(digitalRead(12)==HIGH){
      break;
    }
  }  
  tone(3,1500);
  delay(99999999); 
}