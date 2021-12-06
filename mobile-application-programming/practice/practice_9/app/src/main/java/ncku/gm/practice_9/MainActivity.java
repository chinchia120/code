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
import android.widget.EditText;
import android.widget.TextView;

import com.google.android.gms.maps.CameraUpdate;
import com.google.android.gms.maps.CameraUpdateFactory;
import com.google.android.gms.maps.GoogleMap;
import com.google.android.gms.maps.OnMapReadyCallback;
import com.google.android.gms.maps.SupportMapFragment;
import com.google.android.gms.maps.model.LatLng;
import com.google.android.gms.maps.model.MarkerOptions;

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
        lat_now = location.getLatitude();
        lon_now = location.getLongitude();

        gm.animateCamera(CameraUpdateFactory.newLatLng(new LatLng(lat_now,lon_now)));
        //gm.animateCamera(CameraUpdateFactory.newLatLngZoom(new LatLng(lat_now,lon_now),8f));

        Geocoder geo = new Geocoder(this, Locale.getDefault());
        try {
            String str_now = " ",str_enter = " ";
            List<Address> list_now = geo.getFromLocation(lat_now,lon_now,1);
            List<Address> list_enter = geo.getFromLocationName(((EditText)findViewById(R.id.edt_enter_place)).getText().toString(),1);

            if(list_now!=null && list_now.size()>0){
                Address address_now = list_now.get(0);
                for(int i=0;i<=address_now.getMaxAddressLineIndex();i++){
                    str_now+=address_now.getAddressLine(i);
                }
            }

            if(list_enter!=null && list_enter.size()>0){
                Address address_enter = list_enter.get(0);
                lat_enter = address_enter.getLatitude();
                lon_enter = address_enter.getLongitude();
                str_enter+=lat_enter+"\n"+lon_enter;
            }
            dis = getDistance(lat_now,lon_now,lat_enter,lon_enter)/1000;
            ((TextView)findViewById(R.id.txv_show)).setText(str_now+"\n"+str_enter+"\n"+dis);
            //String.format("距離 : %.02f m",dis)

        }catch (Exception e){
            ((TextView)findViewById(R.id.txv_show)).setText("無法辨識");
        }
        //gm.moveCamera(CameraUpdateFactory.zoomTo(2));
        gm.addMarker(new MarkerOptions().position(new LatLng(lat_now,lon_now)).title("目前位置"));
        gm.addMarker(new MarkerOptions().position(new LatLng(lat_enter,lon_enter)).title("輸入位置"));
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
}