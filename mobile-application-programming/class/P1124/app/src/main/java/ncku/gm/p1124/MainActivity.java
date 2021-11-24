package ncku.gm.p1124;

import androidx.appcompat.app.AppCompatActivity;
import androidx.constraintlayout.widget.ConstraintLayout;

import android.hardware.Sensor;
import android.hardware.SensorEvent;
import android.hardware.SensorEventListener;
import android.hardware.SensorManager;
import android.os.Bundle;
import android.widget.ImageView;
import android.widget.TextView;
import android.widget.Toast;

import java.util.function.ToLongFunction;

public class MainActivity extends AppCompatActivity
        implements SensorEventListener {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        SensorManager snm = (SensorManager) getSystemService(SENSOR_SERVICE);
        Sensor sn_acc = snm.getDefaultSensor(Sensor.TYPE_ACCELEROMETER);
        Sensor sn_mas = snm.getDefaultSensor(Sensor.TYPE_MAGNETIC_FIELD);
        snm.registerListener(this,sn_acc,3);
        snm.registerListener(this,sn_mas,3);
        //snm.registerListener(this,sn_acc,SensorManager.SENSOR_DELAY_NORMAL);
        //Sensor sr = ((SensorManager)getSystemService(SENSOR_SERVICE)).getDefaultSensor(Sensor.TYPE_ACCELEROMETER);
    }

    float [] value_acc = new float[3];
    float [] value_mag = new float[3];
    @Override
    public void onSensorChanged(SensorEvent sensorEvent) {
        ImageView imv_bub = findViewById(R.id.imv_bubble);
        if(sensorEvent.sensor.getType()==Sensor.TYPE_ACCELEROMETER){
            value_acc[0]=sensorEvent.values[0];
            value_acc[1]=sensorEvent.values[1];
            value_acc[2]=sensorEvent.values[2];

            //ConstraintLayout.LayoutParams pr = (ConstraintLayout.LayoutParams) imv_bub.getLayoutParams();
            //pr.horizontalBias = (float) (value_acc[0]*(1/19.6)+0.5);
            //pr.verticalBias = (float) (value_acc[1]*(-1/19.6)+0.5);
            //imv_bub.setLayoutParams(pr);

            float sum = Math.abs(value_acc[0])+Math.abs(value_acc[1])+Math.abs(value_acc[2]);
            if(sum>32){
                Toast.makeText(this,sum+"",Toast.LENGTH_LONG).show();
            }

        }else if(sensorEvent.sensor.getType()==Sensor.TYPE_MAGNETIC_FIELD){
            value_mag[0]=sensorEvent.values[0];
            value_mag[1]=sensorEvent.values[1];
            value_mag[2]=sensorEvent.values[2];

        }
        //((TextView)findViewById(R.id.txv_show)).setText(String.format("X_a : %.8f\nY_a : %.8f\nZ_a : %.8f\nX_m :%.8f\nY_m : %.8f\nZ_m : %.8f",value_acc[0],value_acc[1],value_acc[2],value_mag[0],value_mag[1],value_mag[2]));

        float[] rotation = new float[9];
        float[] degree = new float[3];
        SensorManager.getRotationMatrix(rotation,null,value_acc,value_mag);
        SensorManager.getOrientation(rotation,degree);
        float angle = (float) Math.toDegrees(degree[0]);

        ((TextView)findViewById(R.id.txv_show)).setText(angle+"");
        imv_bub.setRotation(-angle);

    }

    @Override
    public void onAccuracyChanged(Sensor sensor, int i) {

    }
}