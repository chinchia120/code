package ncku.gm.final_project;

import androidx.appcompat.app.AppCompatActivity;

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
import android.widget.Button;
import android.widget.EditText;
import android.widget.TextView;

import com.google.android.material.snackbar.Snackbar;

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

public class MainActivity_show_data_detail extends AppCompatActivity implements View.OnClickListener {

    SQLiteDatabase db;
    UserInformation userInformation = new UserInformation();
    LocationData locationData = new LocationData();
    Cursor cus;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main_show_data_detail);

        ((Button)findViewById(R.id.btn_back_data_detail)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_chat)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_sure)).setOnClickListener(this);

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
    }

    @Override
    public void onClick(View view) {
        if (view.getId()==R.id.btn_back_data_detail){
            finish();
        }else if(view.getId()==R.id.btn_chat){

        }else if(view.getId()==R.id.btn_sure){
            if(cus.getString(1).matches(userInformation.getUser_name())){
                Snackbar.make(findViewById(R.id.root_show_data_detail),"這是自己發起的共乘喔",Snackbar.LENGTH_LONG).show();
            }else {
                startActivity(new Intent(this,MainActivity_together.class));
            }
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
}