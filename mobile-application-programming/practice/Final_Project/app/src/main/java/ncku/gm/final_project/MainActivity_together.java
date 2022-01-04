package ncku.gm.final_project;

import androidx.annotation.NonNull;
import androidx.appcompat.app.AppCompatActivity;
import androidx.core.app.ActivityCompat;

import android.Manifest;
import android.content.Context;
import android.content.pm.PackageManager;
import android.location.Location;
import android.location.LocationManager;
import android.os.Bundle;

import com.google.android.gms.maps.CameraUpdateFactory;
import com.google.android.gms.maps.GoogleMap;
import com.google.android.gms.maps.LocationSource;
import com.google.android.gms.maps.OnMapReadyCallback;
import com.google.android.gms.maps.SupportMapFragment;
import com.google.android.gms.maps.model.LatLng;
import com.google.android.gms.maps.model.MarkerOptions;
import com.google.android.material.snackbar.Snackbar;

public class MainActivity_together extends AppCompatActivity implements LocationSource.OnLocationChangedListener , OnMapReadyCallback {

    GoogleMap gm =null;
    LocationData locationData = new LocationData();

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main_together);

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
        ((LocationManager)getSystemService(Context.LOCATION_SERVICE)).requestLocationUpdates(LocationManager.NETWORK_PROVIDER,1000,5,this::onLocationChanged);
        ((LocationManager)getSystemService(Context.LOCATION_SERVICE)).requestLocationUpdates(LocationManager.GPS_PROVIDER,1000,5,this::onLocationChanged);
        ((SupportMapFragment)getSupportFragmentManager().findFragmentById(R.id.map)).getMapAsync(this);
    }

    @Override
    public void onLocationChanged(@NonNull Location location) {
        if(location!=null){
            gm.animateCamera(CameraUpdateFactory.newLatLngZoom(new LatLng(location.getLatitude(),location.getLongitude()),7f));
            gm.clear();
            gm.addMarker(new MarkerOptions().position(new LatLng(location.getLatitude(),location.getLongitude())).title("目前位置"));
            gm.addMarker(new MarkerOptions().position(new LatLng(locationData.getLat_start(),locationData.getLon_start())).title("出發地"));
            gm.addMarker(new MarkerOptions().position(new LatLng(locationData.getLat_end(),locationData.getLon_end())).title("目的地"));
        }
    }

    @Override
    public void onMapReady(@NonNull GoogleMap googleMap) {
        gm = googleMap;
        gm.setMapType(GoogleMap.MAP_TYPE_NORMAL);

        //Snackbar.make(findViewById(R.id.root_together), locationData.getLat_start() + "\t" + locationData.getLon_start() + "\t" + locationData.getLat_end() + "\t" + locationData.getLon_end(), Snackbar.LENGTH_LONG).show();

    }
}