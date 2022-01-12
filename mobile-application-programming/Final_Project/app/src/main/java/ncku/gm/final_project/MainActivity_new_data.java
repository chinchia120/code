package ncku.gm.final_project;

import androidx.appcompat.app.AppCompatActivity;

import android.app.DatePickerDialog;
import android.app.TimePickerDialog;
import android.content.ContentValues;
import android.content.Context;
import android.content.Intent;
import android.database.Cursor;
import android.database.sqlite.SQLiteDatabase;
import android.location.Address;
import android.location.Geocoder;
import android.location.Location;
import android.os.Bundle;
import android.view.View;
import android.widget.DatePicker;
import android.widget.EditText;
import android.widget.ImageView;
import android.widget.TimePicker;

import com.google.android.material.snackbar.Snackbar;

import java.util.Calendar;
import java.util.List;
import java.util.Locale;
import androidx.annotation.NonNull;
import androidx.appcompat.app.AppCompatActivity;
import androidx.core.app.ActivityCompat;

import android.os.Bundle;
import android.content.pm.PackageManager;
import android.location.Location;
import android.location.LocationListener;
import android.location.LocationManager;
import android.os.Bundle;
import android.os.Handler;
import android.os.Looper;
import android.os.Message;
//import android.support.v4.app.ActivityCompat;
//import android.support.v7.app.AppCompatActivity;
import android.view.View;
import android.widget.Button;
import android.widget.EditText;
import android.widget.TextView;
import android.widget.Toast;

import java.util.HashMap;
import java.util.Random;
import java.util.jar.Manifest;
public class MainActivity_new_data extends AppCompatActivity implements View.OnClickListener , DatePickerDialog.OnDateSetListener , TimePickerDialog.OnTimeSetListener {

    SQLiteDatabase db,db_location;
    UserInformation userInformation = new UserInformation();
    LocationData locationData = new LocationData();

    ClientThread mClientThread;
    private Handler mInputHandler;
    private String host = "140.116.47.94";
    private int port = 7070;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main_new_data);

        db = openOrCreateDatabase("Test_DB", Context.MODE_PRIVATE,null);
        db.execSQL("CREATE TABLE IF NOT EXISTS table_location (_id INTEGER PRIMARY KEY AUTOINCREMENT,name VARCHAR(32),end_place VARCHAR(32),start VARCHAR(32),time VARCHAR(32))");


        ((ImageView)findViewById(R.id.imv_back_user_data)).setOnClickListener(this);
        ((ImageView)findViewById(R.id.imv_finish_new)).setOnClickListener(this);
        ((EditText)findViewById(R.id.edt_time)).setOnClickListener(this);

        mInputHandler = new Handler() {
            @Override
            public void handleMessage(Message msg) {
                if (msg.what == 0) {
                    String[] tmp = msg.obj.toString().split(",");
                    if(tmp[0].matches("location_data")){
                        Cursor cus = db.rawQuery("SELECT * FROM table_location",null);
                        int flag=0;
                        cus.moveToFirst();
                        for(int i=0;i<cus.getCount();i++){
                            if(cus.getString(1).matches(tmp[1]) && cus.getString(2).matches(tmp[2]) && cus.getString(3).matches(tmp[3]) && cus.getString(4).matches(tmp[4])){
                                flag=1;
                                break;
                            }else {
                                cus.moveToNext();
                            }
                        }
                        if (flag==0){
                            ContentValues cv = new ContentValues(4);
                            cv.put("name",tmp[1]);
                            cv.put("end_place",tmp[2]);
                            cv.put("start",tmp[3]);
                            cv.put("time",tmp[4]);
                            db.insert("table_location",null,cv);
                        }
                    }
                }
            }
        };

        mClientThread = new ClientThread(mInputHandler, host, port);
        new Thread(mClientThread).start();
    }



    @Override
    public void onClick(View view) {
        if(view.getId()==R.id.imv_back_user_data){
            finish();
        }else if(view.getId()==R.id.imv_finish_new){
            if(((EditText)findViewById(R.id.edt_start)).getText().toString().matches("") || ((EditText)findViewById(R.id.edt_end)).getText().toString().matches("") || ((EditText)findViewById(R.id.edt_time)).getText().toString().matches("")){
                Snackbar.make(findViewById(R.id.root_new_data),"資訊不完整",Snackbar.LENGTH_LONG).show();
            }else{
                ContentValues cv = new ContentValues(4);
                cv.put("name",UserInformation.getUser_name());
                cv.put("end_place",((EditText)findViewById(R.id.edt_end)).getText().toString());
                cv.put("start",((EditText)findViewById(R.id.edt_start)).getText().toString());
                cv.put("time",((EditText)findViewById(R.id.edt_time)).getText().toString());
                db.insert("table_location",null,cv);

                Cursor cus = db.rawQuery("SELECT * FROM table_location",null);
                cus.moveToFirst();
                for (int i=0;i<cus.getCount();i++){
                    Message msg = new Message();
                    msg.what = 1;
                    msg.obj = "location_data"+","+
                            cus.getString(1)+","+
                            cus.getString(2)+","+
                            cus.getString(3)+","+
                            cus.getString(4);
                    mClientThread.mOutputHandler.sendMessage(msg);
                    cus.moveToNext();
                }
                finish();
            }
        }else if(view.getId()==R.id.edt_time){
            ((EditText)findViewById(R.id.edt_time)).setText("");

            Calendar calendar = Calendar.getInstance();
            int nowyear = calendar.get(Calendar.YEAR);
            int nowmon = calendar.get(Calendar.MONTH);
            int nowdate = calendar.get(Calendar.DATE);


            new DatePickerDialog(this,this,nowyear,nowmon,nowdate).show();
        }
    }


    @Override
    public void onDateSet(DatePicker datePicker, int i, int i1, int i2) {
        Calendar calendar = Calendar.getInstance();
        int nowhour = calendar.get(Calendar.HOUR_OF_DAY);
        int nowmin = calendar.get(Calendar.MINUTE);
        ((EditText)findViewById(R.id.edt_time)).setText(String.format("%02d/%02d",(i1+1),i2));
        new TimePickerDialog(this,this,nowhour,nowmin,true).show();
    }

    @Override
    public void onTimeSet(TimePicker timePicker, int i, int i1) {
        ((EditText)findViewById(R.id.edt_time)).append(String.format(" %02d:%02d",i,i1));
    }
}