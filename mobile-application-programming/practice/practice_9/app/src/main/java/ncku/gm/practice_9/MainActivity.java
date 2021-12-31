package ncku.gm.practice_9;

import androidx.annotation.NonNull;
import androidx.appcompat.app.AppCompatActivity;
import androidx.core.app.ActivityCompat;

import android.Manifest;
import android.app.SearchManager;
import android.content.Context;
import android.content.Intent;
import android.content.pm.PackageManager;
import android.hardware.Sensor;
import android.hardware.SensorEvent;
import android.hardware.SensorEventListener;
import android.hardware.SensorManager;
import android.location.Address;
import android.location.Geocoder;
import android.location.Location;
import android.location.LocationListener;
import android.location.LocationManager;
import android.media.Image;
import android.os.Bundle;
import android.view.View;
import android.widget.Button;
import android.widget.EditText;
import android.widget.ListView;
import android.widget.TextView;
import android.widget.Toast;

import com.google.android.gms.maps.CameraUpdateFactory;
import com.google.android.gms.maps.GoogleMap;
import com.google.android.gms.maps.OnMapReadyCallback;
import com.google.android.gms.maps.SupportMapFragment;
import com.google.android.gms.maps.model.BitmapDescriptor;
import com.google.android.gms.maps.model.BitmapDescriptorFactory;
import com.google.android.gms.maps.model.GroundOverlayOptions;
import com.google.android.gms.maps.model.LatLng;
import com.google.android.gms.maps.model.MarkerOptions;
import com.google.android.material.snackbar.Snackbar;

import java.io.IOException;
import java.util.List;
import java.util.Locale;

public class MainActivity extends AppCompatActivity
        implements LocationListener , OnMapReadyCallback , View.OnClickListener , SensorEventListener {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        if (ActivityCompat.checkSelfPermission(this, Manifest.permission.ACCESS_FINE_LOCATION) != PackageManager.PERMISSION_GRANTED && ActivityCompat.checkSelfPermission(this, Manifest.permission.ACCESS_COARSE_LOCATION) != PackageManager.PERMISSION_GRANTED) {
            // TODO: Consider calling
            //    ActivityCompat#requestPermissions
            ActivityCompat.requestPermissions(this,new String[]{Manifest.permission.ACCESS_FINE_LOCATION,Manifest.permission.ACCESS_COARSE_LOCATION},000);
            // here to request the missing permissions, and then overriding
            //   public void onRequestPermissionsResult(int requestCode, String[] permissions,
            //                                          int[] grantResults)
            // to handle the case where the user grants the permission. See the documentation
            // for ActivityCompat#requestPermissions for more details.
            return;
        }

        ((LocationManager) getSystemService(Context.LOCATION_SERVICE)).requestLocationUpdates("network", 1000, 5, this);
        ((SupportMapFragment)getSupportFragmentManager().findFragmentById(R.id.map)).getMapAsync(this);
        ((Button)findViewById(R.id.btn_search)).setOnClickListener(this);

        SensorManager snm = (SensorManager) getSystemService(SENSOR_SERVICE);
        Sensor sn_acc = snm.getDefaultSensor(Sensor.TYPE_ACCELEROMETER);
        Sensor sn_mas = snm.getDefaultSensor(Sensor.TYPE_MAGNETIC_FIELD);
        snm.registerListener(this,sn_acc,0);
        snm.registerListener(this,sn_mas,0);
    }

    double lat_now,lon_now,lat_enter,lon_enter,dis;
    GoogleMap gm = null;

    @Override
    public void onLocationChanged(@NonNull Location location) {
        if(location!=null){
            lat_now = location.getLatitude();
            lon_now = location.getLongitude();
            Geocoder geo = new Geocoder(this, Locale.TRADITIONAL_CHINESE);
            try {
                List<Address> list_now = geo.getFromLocation(lat_now,lon_now,1);
                String str_now = list_now.get(0).getAddressLine(0);
                ((TextView)findViewById(R.id.txv_show_now)).setText(String.format("目前位置 : %s",str_now));
            } catch (IOException e) {
                e.printStackTrace();
            }

            a = getAngle(lat_now,lon_now,lat_enter,lon_enter);
            if(gm!=null){
                if(dis<20){
                    gm.animateCamera(CameraUpdateFactory.newLatLngZoom(new LatLng(lat_now,lon_now),12f));
                    //gm.addMarker(new MarkerOptions().position(new LatLng(lat_now,lon_now)).anchor(0.5f,0.5f).rotation((float)a+180f).title("輸入位置"));
                    gm.addGroundOverlay(new GroundOverlayOptions().image(BitmapDescriptorFactory.fromResource(R.drawable.arrow)).position(new LatLng(lat_now,lon_now),860f,650f));
                }else if(dis<50){
                    gm.animateCamera(CameraUpdateFactory.newLatLngZoom(new LatLng(lat_now,lon_now),10f));
                    //gm.addMarker(new MarkerOptions().position(new LatLng(lat_now,lon_now)).anchor(0.5f,0.5f).rotation((float)a+180f).title("輸入位置"));
                    gm.addGroundOverlay(new GroundOverlayOptions().image(BitmapDescriptorFactory.fromResource(R.drawable.arrow)).position(new LatLng(lat_now,lon_now),5000f,3000f));
                }else{
                    gm.animateCamera(CameraUpdateFactory.newLatLngZoom(new LatLng(lat_now,lon_now),7f));
                    //gm.addMarker(new MarkerOptions().position(new LatLng(lat_now,lon_now)).anchor(0.5f,0.5f).rotation((float)a+180f).title("輸入位置"));
                    gm.addGroundOverlay(new GroundOverlayOptions().image(BitmapDescriptorFactory.fromResource(R.drawable.arrow)).position(new LatLng(lat_now,lon_now),20000f,18000f));
                }
            }
        }else{
            ((TextView)findViewById(R.id.txv_show_now)).setText("無法取得定位資訊");
        }
    }

    @Override
    public void onMapReady(@NonNull GoogleMap googleMap) {
        gm = googleMap;
        gm.setMapType(GoogleMap.MAP_TYPE_NORMAL);
        Snackbar.make(findViewById(R.id.root),"開始定位",Snackbar.LENGTH_SHORT).show();
    }

    public double getDistance(double lat_now,double lon_now,double lat_enter, double lon_enter){
        float[] result = new float[1];
        Location.distanceBetween(lat_now,lon_now,lat_enter,lon_enter,result);
        return result[0];
    }

    public void onEditTextChange(View view){
        String str_direction = "";
        Geocoder geo = new Geocoder(this, Locale.TRADITIONAL_CHINESE);
        try {
            String str_enter = null;
            if(((EditText)findViewById(R.id.edt_enter_place)).getText().toString()!=null){
                List<Address> list_enter = geo.getFromLocationName(((EditText)findViewById(R.id.edt_enter_place)).getText().toString(),1);
                lat_enter = list_enter.get(0).getLatitude();
                lon_enter = list_enter.get(0).getLongitude();
                if(lat_enter==0 || lon_enter==0){
                    Snackbar.make(findViewById(R.id.root),"找不到適合目的",Snackbar.LENGTH_SHORT).show();
                    return;
                }
                List<Address> list_enter_place = geo.getFromLocation(lat_enter,lon_enter,1);
                str_enter = list_enter_place.get(0).getAddressLine(0);
                dis = getDistance(lat_now,lon_now,lat_enter,lon_enter)/1000;
                if(lon_enter>lon_now){
                    str_direction+="東";
                }else{
                    str_direction+="西";
                }
                if(lat_enter>lat_now){
                    str_direction+="北";
                }else{
                    str_direction+="南";
                }
            }else{
                str_enter = "無法辨識輸入地點";
                dis = 0f;
            }
            //((TextView)findViewById(R.id.txv_show_enter)).setText(String.format("%f\n%f",lat_enter,lon_enter));
            ((TextView)findViewById(R.id.txv_show_enter)).setText(String.format("輸入位置 : %s",str_enter));
            ((TextView)findViewById(R.id.txv_show_dis)).setText(String.format("距離 : %.1f km , %s方",dis,str_direction));
            Snackbar.make(findViewById(R.id.root),"輸入成功",Snackbar.LENGTH_SHORT).show();
            if(gm!=null){
                gm.clear();
                if(dis<20){
                    gm.animateCamera(CameraUpdateFactory.newLatLngZoom(new LatLng(lat_now,lon_now),12f));
                    //gm.addMarker(new MarkerOptions().position(new LatLng(lat_now,lon_now)).anchor(0.5f,0.5f).rotation((float)a+180f).title("輸入位置"));
                    gm.addGroundOverlay(new GroundOverlayOptions().image(BitmapDescriptorFactory.fromResource(R.drawable.arrow)).position(new LatLng(lat_now,lon_now),860f,650f));
                }else if(dis<50){
                    gm.animateCamera(CameraUpdateFactory.newLatLngZoom(new LatLng(lat_now,lon_now),10f));
                    //gm.addMarker(new MarkerOptions().position(new LatLng(lat_now,lon_now)).anchor(0.5f,0.5f).rotation((float)a+180f).title("輸入位置"));
                    gm.addGroundOverlay(new GroundOverlayOptions().image(BitmapDescriptorFactory.fromResource(R.drawable.arrow)).position(new LatLng(lat_now,lon_now),5000f,3000f));
                }else{
                    gm.animateCamera(CameraUpdateFactory.newLatLngZoom(new LatLng(lat_now,lon_now),7f));
                    //gm.addMarker(new MarkerOptions().position(new LatLng(lat_now,lon_now)).anchor(0.5f,0.5f).rotation((float)a+180f).title("輸入位置"));
                    gm.addGroundOverlay(new GroundOverlayOptions().image(BitmapDescriptorFactory.fromResource(R.drawable.arrow)).position(new LatLng(lat_now,lon_now),20000f,18000f));
                }
                gm.addMarker(new MarkerOptions().position(new LatLng(lat_enter,lon_enter)).title("輸入位置"));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Override
    public void onClick(View view) {
        startActivity(new Intent(Intent.ACTION_WEB_SEARCH).putExtra(SearchManager.QUERY,"附近景點"));
    }

    float [] value_acc = new float[3];
    float [] value_mag = new float[3];
    float angle,angle_enter;
    double a;

    @Override
    public void onSensorChanged(SensorEvent sensorEvent) {
        if(sensorEvent.sensor.getType()==Sensor.TYPE_ACCELEROMETER){
            for(int i=0;i<3;i++){
                value_acc[i]=sensorEvent.values[i];
            }
        }else if(sensorEvent.sensor.getType()==Sensor.TYPE_MAGNETIC_FIELD){
            for(int i=0;i<3;i++){
                value_mag[i]=sensorEvent.values[i];
            }
        }
        float[] rotation = new float[9];
        float[] degree = new float[3];
        SensorManager.getRotationMatrix(rotation,null,value_acc,value_mag);
        SensorManager.getOrientation(rotation,degree);
        //angle = (float) Math.toDegrees(degree[0]);
        //angle_enter = (float) Math.atan(lon_enter/lat_enter);
    }

    @Override
    public void onAccuracyChanged(Sensor sensor, int i) {

    }

    public double getAngle(double lat_a, double lng_a, double lat_b, double lng_b) {

        double d = 0;
        lat_a=lat_a*Math.PI/180;
        lng_a=lng_a*Math.PI/180;
        lat_b=lat_b*Math.PI/180;
        lng_b=lng_b*Math.PI/180;

        d=Math.sin(lat_a)*Math.sin(lat_b)+Math.cos(lat_a)*Math.cos(lat_b)*Math.cos(lng_b-lng_a);
        d=Math.sqrt(1-d*d);
        d=Math.cos(lat_b)*Math.sin(lng_b-lng_a)/d;
        d=Math.asin(d)*180/Math.PI;

        //d = Math.round(d*10000);
        return d;
    }
}