package ncku.gm.practice_8;

import androidx.appcompat.app.AlertDialog;
import androidx.appcompat.app.AppCompatActivity;
import androidx.constraintlayout.widget.ConstraintLayout;

import android.content.Context;
import android.content.DialogInterface;
import android.content.Intent;
import android.hardware.Sensor;
import android.hardware.SensorEvent;
import android.hardware.SensorEventListener;
import android.hardware.SensorManager;
import android.os.Bundle;
import android.os.CountDownTimer;
import android.os.Vibrator;
import android.view.View;
import android.widget.Button;
import android.widget.ImageView;
import android.widget.TextView;

import com.google.android.material.snackbar.Snackbar;

import java.util.Random;

public class MainActivity2 extends AppCompatActivity implements SensorEventListener , View.OnClickListener , DialogInterface.OnClickListener {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main2);

        snm = (SensorManager) getSystemService(SENSOR_SERVICE);
        Sensor sn_acc = snm.getDefaultSensor(Sensor.TYPE_ACCELEROMETER);
        snm.registerListener(this,sn_acc,1);

        ((Button)findViewById(R.id.btn_reset)).setOnClickListener(this);

    }
    SensorManager snm;
    float x=(float)0.50,y=(float)0.20,x_now,y_now;
    float [] value_acc_eat = new float[3];
    int cnt_reset=0,cnt_eat=0,flag=0,flag_next=2;
    Vibrator vb;

    @Override
    public void onSensorChanged(SensorEvent sensorEvent) {
        for(int i=0;i<3;i++){
            value_acc_eat[i] = sensorEvent.values[i];
        }
        if(value_acc_eat[0]>9.8 || value_acc_eat[0]<-9.8 || value_acc_eat[1]>9.8 || value_acc_eat[1]<-9.8){
            ((ImageView)findViewById(R.id.imv_eat)).setVisibility(View.GONE);
        }else{
            ((ImageView)findViewById(R.id.imv_eat)).setVisibility(View.VISIBLE);
        }
        x_now = (float) Math.round((value_acc_eat[0]*(-1/19.6)+0.5)*100)/100;
        y_now = (float) Math.round((value_acc_eat[1]*(1/19.6)+0.5)*100)/100;
        ConstraintLayout.LayoutParams pr = (ConstraintLayout.LayoutParams)((ImageView)findViewById(R.id.imv_eat)).getLayoutParams();
        pr.horizontalBias = x_now;
        pr.verticalBias = y_now;
        ((ImageView)findViewById(R.id.imv_eat)).setLayoutParams(pr);

        //((TextView)findViewById(R.id.txv_show)).setText(String.format("%.2f %.2f\n%.2f %.2f",x,y,x_now,y_now));

        if(x_now==x && y_now==y){
            set_random_number();
            if(flag==0){
                flag=1;
                Snackbar.make(findViewById(R.id.root),"開始遊戲",Snackbar.LENGTH_SHORT).show();
            }else if(flag==1){
                cnt_eat++;
                if(cnt_eat==5){
                    flag_next=1;
                    AlertDialog.Builder bud = new AlertDialog.Builder(this);
                    bud.setTitle("通關成功");
                    bud.setMessage("是否前往下一關");
                    bud.setPositiveButton("確定",this);
                    bud.setNegativeButton("取消",this);
                    bud.setCancelable(false);
                    bud.show();
                }
            }
        }
        ((TextView)findViewById(R.id.txv_show)).setText(String.format("目前累積吃掉%d次\n目前累積重置%d次",cnt_eat,cnt_reset));
    }

    @Override
    public void onAccuracyChanged(Sensor sensor, int i) {

    }

    public void set_random_number(){
        x = (float) (new Random().nextInt(101))/100;
        y = (float) (new Random().nextInt(101))/100;

        ConstraintLayout.LayoutParams pr = (ConstraintLayout.LayoutParams)((ImageView)findViewById(R.id.imv_star)).getLayoutParams();
        pr.horizontalBias = x;
        pr.verticalBias = y ;
        ((ImageView)findViewById(R.id.imv_star)).setLayoutParams(pr);
    }

    @Override
    public void onClick(View view) {
        cnt_reset++;
        if(cnt_reset==3){
            vb =(Vibrator) getSystemService(Context.VIBRATOR_SERVICE);
            vb.vibrate(1000);
            flag_next=0;
            AlertDialog.Builder bud = new AlertDialog.Builder(this);
            bud.setTitle("通關失敗");
            bud.setMessage("是否重新開始");
            bud.setPositiveButton("確定",this);
            bud.setNegativeButton("取消",this);
            bud.setCancelable(false);
            bud.show();
        }
        set_random_number();
    }

    @Override
    public void onClick(DialogInterface dialogInterface, int i) {
        if(flag_next==1){
            if(i==-1){
                snm.unregisterListener(this);
                startActivity(new Intent(this,MainActivity3.class));
            }else if(i==-2){
                snm.unregisterListener(this);
                startActivity(new Intent(this,MainActivity4.class));
            }
        }else if(flag_next==0){
            if(i==-1){
                cnt_eat=0;
                cnt_reset=0;
                x=(float)0.50;
                y=(float)0.20;
                ((TextView)findViewById(R.id.txv_show)).setText(String.format("目前累積吃掉%d次\n目前累積重置%d次",cnt_eat,cnt_reset));
                ConstraintLayout.LayoutParams pr = (ConstraintLayout.LayoutParams)((ImageView)findViewById(R.id.imv_star)).getLayoutParams();
                pr.horizontalBias = x;
                pr.verticalBias = y ;
                ((ImageView)findViewById(R.id.imv_star)).setLayoutParams(pr);
            }else if(i==-2){
                snm.unregisterListener(this);
                startActivity(new Intent(this,MainActivity4.class));
            }
        }
    }
}