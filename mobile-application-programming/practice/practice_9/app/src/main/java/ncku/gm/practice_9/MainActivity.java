package ncku.gm.practice_9;

import androidx.annotation.NonNull;
import androidx.appcompat.app.AppCompatActivity;
import androidx.core.app.ActivityCompat;

import android.Manifest;
import android.content.Context;
import android.content.pm.PackageManager;
import android.location.Address;
import android.location.Geocoder;
import android.location.Location;
import android.location.LocationListener;
import android.location.LocationManager;
import android.os.Bundle;
import android.view.View;
import android.widget.EditText;
import android.widget.TextView;

import com.google.android.gms.maps.CameraUpdateFactory;
import com.google.android.gms.maps.GoogleMap;
import com.google.android.gms.maps.OnMapReadyCallback;
import com.google.android.gms.maps.SupportMapFragment;
import com.google.android.gms.maps.model.LatLng;
import com.google.android.gms.maps.model.MarkerOptions;

import java.io.IOException;
import java.util.List;
import java.util.Locale;

public class MainActivity extends AppCompatActivity
        implements LocationListener , OnMapReadyCallback {

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
        ((LocationManager) getSystemService(Context.LOCATION_SERVICE)).requestLocationUpdates("network", 5000, 5, this);

        ((SupportMapFragment)getSupportFragmentManager().findFragmentById(R.id.map)).getMapAsync(this);
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
            if(gm!=null){
                gm.animateCamera(CameraUpdateFactory.newLatLngZoom(new LatLng(lat_now,lon_now),7f));
                gm.addMarker(new MarkerOptions().position(new LatLng(lat_now,lon_now)).title("目前位置"));
            }
        }else{
            ((TextView)findViewById(R.id.txv_show_now)).setText("無法取得定位資訊");
        }
    }

    @Override
    public void onMapReady(@NonNull GoogleMap googleMap) {
        gm = googleMap;
        gm.setMapType(GoogleMap.MAP_TYPE_NORMAL);
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
            ((TextView)findViewById(R.id.txv_show_enter)).setText(String.format("輸入位置 : %s",str_enter));
            ((TextView)findViewById(R.id.txv_show_dis)).setText(String.format("距離 : %.1f km , %s方",dis,str_direction));
            if(gm!=null){
                gm.addMarker(new MarkerOptions().position(new LatLng(lat_enter,lon_enter)).title("輸入位置"));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}