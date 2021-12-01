package ncku.gm.p1201;

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
import android.widget.Button;
import android.widget.EditText;
import android.widget.TextView;

import com.google.android.gms.maps.CameraUpdateFactory;
import com.google.android.gms.maps.GoogleMap;
import com.google.android.gms.maps.OnMapReadyCallback;
import com.google.android.gms.maps.SupportMapFragment;
import com.google.android.gms.maps.model.LatLng;

import java.util.List;
import java.util.Locale;

public class MainActivity extends AppCompatActivity
        implements LocationListener , View.OnClickListener , OnMapReadyCallback {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        ((Button)findViewById(R.id.btn_address)).setOnClickListener(this);

        LocationManager lm = (LocationManager) getSystemService(Context.LOCATION_SERVICE);
        if (ActivityCompat.checkSelfPermission(this, Manifest.permission.ACCESS_FINE_LOCATION) != PackageManager.PERMISSION_GRANTED && ActivityCompat.checkSelfPermission(this, Manifest.permission.ACCESS_COARSE_LOCATION) != PackageManager.PERMISSION_GRANTED) {
            // TODO: Consider calling
            //    ActivityCompat#requestPermissions
            ActivityCompat.requestPermissions(this,new String[] {Manifest.permission.ACCESS_COARSE_LOCATION,Manifest.permission.ACCESS_FINE_LOCATION},000);
            // here to request the missing permissions, and then overriding
            //   public void onRequestPermissionsResult(int requestCode, String[] permissions,
            //                                          int[] grantResults)
            // to handle the case where the user grants the permission. See the documentation
            // for ActivityCompat#requestPermissions for more details.
            return;
        }
        lm.requestLocationUpdates("gps", 5000, 5, this);

        SupportMapFragment smf = (SupportMapFragment) getSupportFragmentManager().findFragmentById(R.id.map);
        smf.getMapAsync(this);

    }

    double lat=0d,lon=0d,alt=0d;

    @Override
    public void onLocationChanged(@NonNull Location location) {
        lat = location.getLatitude();
        lon = location.getLongitude();
        alt = location.getAltitude();
        ((EditText)findViewById(R.id.edt_lat)).setText(lat+"");
        ((EditText)findViewById(R.id.edt_lon)).setText(lon+"");
        ((EditText)findViewById(R.id.edt_alt)).setText(alt+"");

    }

    @Override
    public void onClick(View view) {

        Geocoder geo = new Geocoder(this, Locale.getDefault());
        try {
            List<Address> la = geo.getFromLocation(lat,lon,1);
            if(la!=null && la.size()>0){
                Address a= la.get(0);
                String s = " ";
                for(int i=0;i<=a.getMaxAddressLineIndex();i++){
                    s+=a.getAddressLine(i);
                }
                ((TextView)findViewById(R.id.txv_show)).setText(s);
            }
        }
        catch (Exception e){

        }
    }

    GoogleMap gm = null;

    @Override
    public void onMapReady(GoogleMap googleMap) {
        gm = googleMap;
        gm.setMapType(GoogleMap.MAP_TYPE_NORMAL);

        //gm.moveCamera(CameraUpdateFactory.zoomTo(18));
        if(gm != null){
            gm.animateCamera(CameraUpdateFactory.newLatLng(new LatLng(lat,lon)));
        }
    }
}