void setup(){
  Serial.begin(9600);
  pinMode(11,OUTPUT);
  pinMode(10,INPUT);
}

void loop(){
  float duration,distance;
  digitalWrite(11,HIGH);
  delayMicroseconds(10);
  digitalWrite(11,LOW);
  duration=pulseIn(10,HIGH); 
  distance=(duration/2)*0.034;
  if(distance<8){
    tone(12,494,200);
  }
}