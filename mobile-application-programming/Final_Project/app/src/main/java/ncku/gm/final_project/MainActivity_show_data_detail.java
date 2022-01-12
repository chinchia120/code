package ncku.gm.final_project;

import androidx.appcompat.app.AlertDialog;
import androidx.appcompat.app.AppCompatActivity;

import android.content.Context;
import android.content.DialogInterface;
import android.content.Intent;
import android.database.Cursor;
import android.database.sqlite.SQLiteDatabase;
import android.hardware.Sensor;
import android.hardware.SensorEvent;
import android.hardware.SensorEventListener;
import android.hardware.SensorManager;
import android.location.Address;
import android.location.Geocoder;
import android.location.Location;
import android.os.Bundle;
import android.os.Handler;
import android.os.Message;
import android.view.View;
import android.widget.ArrayAdapter;
import android.widget.Button;
import android.widget.ListView;
import android.widget.TextView;

import com.google.android.material.snackbar.Snackbar;

import java.util.Calendar;
import java.util.List;
import java.util.Locale;

//import android.support.v4.app.ActivityCompat;
//import android.support.v7.app.AppCompatActivity;

public class MainActivity_show_data_detail extends AppCompatActivity implements View.OnClickListener , SensorEventListener , DialogInterface.OnClickListener {

    SQLiteDatabase db;
    UserInformation userInformation = new UserInformation();
    LocationData locationData = new LocationData();
    Cursor cus;


    ClientThread mClientThread;
    private Handler mInputHandler;
    private String host = "140.116.47.94";
    private int port = 7070;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main_show_data_detail);

        ((Button)findViewById(R.id.btn_back_data_detail)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_chat)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_sure)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_delete)).setOnClickListener(this);

        db = openOrCreateDatabase("Test_DB", Context.MODE_PRIVATE,null);
        db.execSQL("CREATE TABLE IF NOT EXISTS table_location (_id INTEGER PRIMARY KEY AUTOINCREMENT,name VARCHAR(32),end_place VARCHAR(32),start VARCHAR(32),time VARCHAR(32))");

        cus = db.rawQuery("SELECT * FROM table_location",null);
        if(cus.moveToFirst()){
            int i = getIntent().getBundleExtra("Bundle").getInt("Cnt");
            while (i>0){
                cus.moveToNext();
                i--;
            }
        }
        getCoordinate();

        ((TextView)findViewById(R.id.txv_show_name)).setText("發起人 : "+cus.getString(1));
        ((TextView)findViewById(R.id.txv_show_end)).setText(cus.getString(2));
        ((TextView)findViewById(R.id.txv_show_start)).setText("出發處 : "+cus.getString(3));
        ((TextView)findViewById(R.id.txv_show_distance)).setText("距離 : "+LocationData.getDis()+" km");
        ((TextView)findViewById(R.id.txv_show_time)).setText("時間 : "+cus.getString(4));

        //取得系統感測器服務
        SensorManager sm = (SensorManager) getSystemService(SENSOR_SERVICE);
        //透過getDefaultSensor取得裝置預設感測器
        Sensor sr1 = sm.getDefaultSensor(Sensor.TYPE_ACCELEROMETER);
        Sensor sr2 = sm.getDefaultSensor(Sensor.TYPE_MAGNETIC_FIELD);
        //透過registerListener註冊監聽器在感測器上
        sm.registerListener(this,sr1,SensorManager.SENSOR_DELAY_NORMAL);
        sm.registerListener(this,sr2,SensorManager.SENSOR_DELAY_NORMAL);

        if (userInformation.getUser_name().matches(cus.getString(1))){
            ((Button)findViewById(R.id.btn_delete)).setVisibility(View.VISIBLE);
            ((Button)findViewById(R.id.btn_sure)).setVisibility(View.GONE);
        }
        else {
            ((Button)findViewById(R.id.btn_delete)).setVisibility(View.GONE);
            ((Button)findViewById(R.id.btn_sure)).setVisibility(View.VISIBLE);
        }

        mInputHandler = new Handler()
        {
            @Override
            public void handleMessage(Message msg)
            {
                if (msg.what == 0)
                {
                    String[] tmp = msg.obj.toString().split(",");
                }
            }
        };

        mClientThread = new ClientThread(mInputHandler, host, port);
        new Thread(mClientThread).start();
    }

    @Override
    public void onClick(View view) {
        if (view.getId()==R.id.btn_back_data_detail){
            finish();
        }else if(view.getId()==R.id.btn_chat){
            Intent it = new Intent(this,MainActivity_chat.class);
            Bundle bdl = new Bundle();
            bdl.putString("name",cus.getString(1));
            it.putExtra("Bundle",bdl);
            startActivity(it);
        }else if(view.getId()==R.id.btn_sure){
                Message msg = new Message();
                msg.what = 1;
                msg.obj = "together"+","+
                        cus.getString(1)+","+
                        locationData.getLat_start()+","+
                        locationData.getLon_start()+","+
                        locationData.getLat_end()+","+
                        locationData.getLon_end();
                mClientThread.mOutputHandler.sendMessage(msg);

                Message msg_del = new Message();
                msg_del.what = 1;
                msg_del.obj = "delete"+","+
                        cus.getString(1)+","+
                        cus.getString(2)+","+
                        cus.getString(3)+","+
                        cus.getString(4);
                mClientThread.mOutputHandler.sendMessage(msg_del);

                Intent it = new Intent(this,MainActivity_together.class);
                Bundle bdl = new Bundle();
                bdl.putString("id",cus.getString(0));
                it.putExtra("Bundle",bdl);
                startActivity(it);

        }else if (view.getId()==R.id.btn_delete){
            AlertDialog.Builder b = new AlertDialog.Builder(this);
            b.setIcon(android.R.drawable.presence_away);
            b.setTitle("刪除此筆共乘活動");
            b.setMessage("確定刪除?");
            b.setPositiveButton("確定",this);
            b.setNeutralButton("取消",this);
            b.setCancelable(false);
            b.show();
        }
    }

    public double getDistance(double lat_start,double lon_start,double lat_end, double lon_end){
        float[] result = new float[1];
        Location.distanceBetween(lat_start,lon_start,lat_end,lon_end,result);
        return result[0];
    }

    public void getCoordinate(){
        Geocoder geo = new Geocoder(MainActivity_show_data_detail.this, Locale.TRADITIONAL_CHINESE);
        try {
            List<Address> list_start = geo.getFromLocationName(cus.getString(3),1);
            locationData.setLat_start(list_start.get(0).getLatitude());
            locationData.setLon_start(list_start.get(0).getLongitude());

            List<Address> list_end = geo.getFromLocationName(cus.getString(2),1);
            locationData.setLat_end(list_end.get(0).getLatitude());
            locationData.setLon_end(list_end.get(0).getLongitude());

            locationData.setDis(String.format("%.1f",getDistance(locationData.getLat_start(),locationData.getLon_start(),locationData.getLat_end(),locationData.getLon_end())/1000));
        }catch (Exception e){
            e.printStackTrace();
        }
    }


    @Override
    public void onSensorChanged(SensorEvent sensorEvent) {
        Calendar calendar = Calendar.getInstance();

        int nowmon = calendar.get(Calendar.MONTH)+1;
        int nowdate = calendar.get(Calendar.DATE);
        //計算剩下時間
        int starttime = Integer.parseInt(cus.getString(4).substring(6,8)) * 3600
                + Integer.parseInt(cus.getString(4).substring(9,11)) * 60 ;

        int now = calendar.get(Calendar.HOUR_OF_DAY) * 3600
                + calendar.get(Calendar.MINUTE) * 60
                + calendar.get(Calendar.SECOND);

        int lefttime = Math.max(starttime - now,0);

        int startmonth = Integer.parseInt(cus.getString(4).substring(0,2));
        int startdate = Integer.parseInt(cus.getString(4).substring(3,5));

        if (nowmon == startmonth && nowdate == startdate){
            ((TextView)findViewById(R.id.tvlefttime)).setText("距出發時間剩餘 : "
                    + String.format("%02d時%02d分%02d秒", lefttime / 3600, (lefttime % 3600) / 60, lefttime % 60));
        }
        else {
            ((TextView)findViewById(R.id.tvlefttime)).setVisibility(View.GONE);
        }
    }

    @Override
    public void onAccuracyChanged(Sensor sensor, int i) {

    }

    @Override
    public void onClick(DialogInterface dialogInterface, int i) {
        if (i == DialogInterface.BUTTON_POSITIVE){
            Message msg = new Message();
            msg.what = 1;
            msg.obj = "delete"+","+
                    cus.getString(1)+","+
                    cus.getString(2)+","+
                    cus.getString(3)+","+
                    cus.getString(4);
            mClientThread.mOutputHandler.sendMessage(msg);
            finish();
        }
        else if (i == DialogInterface.BUTTON_NEUTRAL){

        }
    }
}