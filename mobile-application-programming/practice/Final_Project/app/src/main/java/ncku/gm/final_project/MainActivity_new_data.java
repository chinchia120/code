package ncku.gm.final_project;

import androidx.appcompat.app.AppCompatActivity;

import android.app.DatePickerDialog;
import android.app.TimePickerDialog;
import android.content.Intent;
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

import java.util.List;
import java.util.Locale;

public class MainActivity_new_data extends AppCompatActivity implements View.OnClickListener , DatePickerDialog.OnDateSetListener , TimePickerDialog.OnTimeSetListener {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main_new_data);

        ((ImageView)findViewById(R.id.imv_back_new)).setOnClickListener(this);
        ((ImageView)findViewById(R.id.imv_finish_new)).setOnClickListener(this);
        ((EditText)findViewById(R.id.edt_time)).setOnClickListener(this);
    }


    @Override
    public void onClick(View view) {
        if(view.getId()==R.id.imv_back_new){
            finish();
        }else if(view.getId()==R.id.imv_finish_new){
            if(((EditText)findViewById(R.id.edt_start)).getText().toString().matches("") || ((EditText)findViewById(R.id.edt_end)).getText().toString().matches("") || ((EditText)findViewById(R.id.edt_time)).getText().toString().matches("")){
                Snackbar.make(findViewById(R.id.root_new_data),"資訊不完整",Snackbar.LENGTH_LONG).show();
            }else{
                double lat_start,lon_start,lat_end,lon_end;
                String str_distance = null;
                Geocoder geo = new Geocoder(this, Locale.TRADITIONAL_CHINESE);
                try {
                    List<Address> list_enter = geo.getFromLocationName(((EditText)findViewById(R.id.edt_start)).getText().toString(),1);
                    lat_start = list_enter.get(0).getLatitude();
                    lon_start = list_enter.get(0).getLongitude();
                    List<Address> list_end = geo.getFromLocationName(((EditText)findViewById(R.id.edt_end)).getText().toString(),1);
                    lat_end = list_end.get(0).getLatitude();
                    lon_end = list_end.get(0).getLongitude();
                    str_distance = String.format("%.1f km",getDistance(lat_start,lon_start,lat_end,lon_end)/1000);
                }catch (Exception e){

                }

                Intent it = new Intent();
                it.putExtra("End",((EditText)findViewById(R.id.edt_end)).getText().toString());
                it.putExtra("Start",((EditText)findViewById(R.id.edt_start)).getText().toString());
                it.putExtra("Time",((EditText)findViewById(R.id.edt_time)).getText().toString());
                it.putExtra("Distance",str_distance);
                setResult(-1,it);
                finish();
            }
        }else if(view.getId()==R.id.edt_time){
            ((EditText)findViewById(R.id.edt_time)).setText("");
            new DatePickerDialog(this,this,2022,0,1).show();
        }
    }

    @Override
    public void onDateSet(DatePicker datePicker, int i, int i1, int i2) {
        ((EditText)findViewById(R.id.edt_time)).setText((i1+1)+"/"+(i2));
        new TimePickerDialog(this,this,12,00,true).show();
    }

    @Override
    public void onTimeSet(TimePicker timePicker, int i, int i1) {
        ((EditText)findViewById(R.id.edt_time)).append(String.format(" %02d:%02d\n",i,i1));
    }

    public double getDistance(double lat_start,double lon_start,double lat_end, double lon_end){
        float[] result = new float[1];
        Location.distanceBetween(lat_start,lon_start,lat_end,lon_end,result);
        return result[0];
    }
}