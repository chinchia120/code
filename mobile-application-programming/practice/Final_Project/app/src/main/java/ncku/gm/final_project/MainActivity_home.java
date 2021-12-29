package ncku.gm.final_project;

import android.os.Bundle;
import android.view.View;
import android.view.Menu;

import com.google.android.material.snackbar.Snackbar;
import com.google.android.material.navigation.NavigationView;

import androidx.navigation.NavController;
import androidx.navigation.Navigation;
import androidx.navigation.ui.AppBarConfiguration;
import androidx.navigation.ui.NavigationUI;
import androidx.drawerlayout.widget.DrawerLayout;
import androidx.appcompat.app.AppCompatActivity;

import ncku.gm.final_project.databinding.ActivityMainBinding;
import ncku.gm.final_project.databinding.ActivityMainHomeBinding;
import android.Manifest;
import android.content.Intent;
import android.content.pm.PackageManager;
import android.location.Location;
import android.location.LocationListener;
import android.location.LocationManager;
import android.os.Bundle;
import android.view.View;
import android.view.Menu;
import android.widget.ArrayAdapter;
import android.widget.Button;

import com.google.android.gms.maps.CameraUpdateFactory;
import com.google.android.gms.maps.GoogleMap;
import com.google.android.gms.maps.OnMapReadyCallback;
import com.google.android.gms.maps.SupportMapFragment;
import com.google.android.gms.maps.model.LatLng;
import com.google.android.gms.maps.model.MarkerOptions;
import com.google.android.material.snackbar.Snackbar;
import com.google.android.material.navigation.NavigationView;

import androidx.annotation.NonNull;
import androidx.annotation.Nullable;
import androidx.core.app.ActivityCompat;
import androidx.navigation.NavController;
import androidx.navigation.Navigation;
import androidx.navigation.ui.AppBarConfiguration;
import androidx.navigation.ui.NavigationUI;
import androidx.drawerlayout.widget.DrawerLayout;
import androidx.appcompat.app.AppCompatActivity;

import java.io.Serializable;
import java.util.ArrayList;

public class MainActivity_home extends AppCompatActivity implements OnMapReadyCallback , LocationListener , View.OnClickListener {

    private AppBarConfiguration mAppBarConfiguration;
    private ActivityMainHomeBinding binding;

    private GoogleMap mMap;
    LatLng myplace;
    double lat;
    double lng;
    ArrayList<String> str_end = new ArrayList<>();
    ArrayList<String> str_start = new ArrayList<>();
    ArrayList<String> str_distance = new ArrayList<>();
    ArrayList<String> str_time = new ArrayList<>();



    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);

//側拉頁面
        binding = ActivityMainHomeBinding.inflate(getLayoutInflater());
        setContentView(binding.getRoot());

        setSupportActionBar(binding.appBarMain.toolbar);
        binding.appBarMain.toolbar.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                Snackbar.make(view, "Replace with your own action", Snackbar.LENGTH_LONG)
                        .setAction("Action", null).show();
            }
        });
        DrawerLayout drawer = binding.drawerLayout;
        //NavigationView navigationView = binding.navView;
        // Passing each menu ID as a set of Ids because each
        // menu should be considered as top level destinations.
        mAppBarConfiguration = new AppBarConfiguration.Builder(
                R.id.nav_home, R.id.nav_gallery, R.id.nav_slideshow)
                .setOpenableLayout(drawer)
                .build();
        NavController navController = Navigation.findNavController(this, R.id.nav_host_fragment_content_main);
        NavigationUI.setupActionBarWithNavController(this, navController, mAppBarConfiguration);
        //NavigationUI.setupWithNavController(navigationView, navController);

        //設定位置
        LocationManager lm = (LocationManager) getSystemService(LOCATION_SERVICE);
        if (ActivityCompat.checkSelfPermission(this, Manifest.permission.ACCESS_FINE_LOCATION) != PackageManager.PERMISSION_GRANTED
                && ActivityCompat.checkSelfPermission(this, Manifest.permission.ACCESS_COARSE_LOCATION) != PackageManager.PERMISSION_GRANTED) {
            // TODO: Consider calling
            //    ActivityCompat#requestPermissions
            ActivityCompat.requestPermissions(this, new String[]{
                    Manifest.permission.ACCESS_COARSE_LOCATION,Manifest.permission.ACCESS_FINE_LOCATION},200);
            // here to request the missing permissions, and then overriding
            //   public void onRequestPermissionsResult(int requestCode, String[] permissions,
            //                                          int[] grantResults)
            // to handle the case where the user grants the permission. See the documentation
            // for ActivityCompat#requestPermissions for more details.
            return;
        }
        lm.requestLocationUpdates("network", 5000, 5, this);

        SupportMapFragment smf = (SupportMapFragment) getSupportFragmentManager().findFragmentById(R.id.map);
        smf.getMapAsync(this);

        ((Button)findViewById(R.id.btn_new_data)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_show_data)).setOnClickListener(this);

        str_end.add("成大醫院");
        str_end.add("台南高鐵站");

        str_start.add("台南火車站");
        str_start.add("自強校區");

        str_time.add("1/1 12:00");
        str_time.add("1/3 17:20");

        str_distance.add("0.9 km");
        str_distance.add("10.4 km");
    }

    @Override
    public boolean onCreateOptionsMenu(Menu menu) {
        // Inflate the menu; this adds items to the action bar if it is present.
        getMenuInflater().inflate(R.menu.main_activity_home, menu);
        return true;
    }

    @Override
    public boolean onSupportNavigateUp() {
        NavController navController = Navigation.findNavController(this, R.id.nav_host_fragment_content_main);
        return NavigationUI.navigateUp(navController, mAppBarConfiguration)
                || super.onSupportNavigateUp();
    }

    @Override
    public void onMapReady(@NonNull GoogleMap googleMap) {
        mMap = googleMap;
        mMap.setMapType(GoogleMap.MAP_TYPE_NORMAL);
    }

    @Override
    public void onLocationChanged(@NonNull Location location) {
        lat = location.getLatitude();
        lng = location.getLongitude();

        if (location != null){

            myplace = new LatLng(lat,lng);

            if (mMap != null){
                mMap.animateCamera(CameraUpdateFactory.newLatLngZoom(new LatLng(lat,lng),15f));
                mMap.addMarker(new MarkerOptions().position(myplace).title("目前位置"));
            }
        }
        else {

        }
    }

    @Override
    public void onClick(View view) {
        if(view.getId()==R.id.btn_new_data){
            startActivityForResult(new Intent(this,MainActivity_new_data.class),002);
        }else if(view.getId()==R.id.btn_show_data){
            Intent it = new Intent(this,MainActivity_show_data.class);
            Bundle bdl = new Bundle();
            bdl.putSerializable("End",(Serializable) str_end);
            bdl.putSerializable("Start",(Serializable) str_start);
            bdl.putSerializable("Time",(Serializable) str_time);
            bdl.putSerializable("Distance",(Serializable) str_distance);
            it.putExtra("Bundle",bdl);
            startActivity(it);
        }
    }

    @Override
    protected void onActivityResult(int requestCode, int resultCode, @Nullable Intent data) {
        super.onActivityResult(requestCode, resultCode, data);

        if(requestCode==002 && resultCode==-1){
            str_end.add(data.getStringExtra("End"));
            str_start.add(data.getStringExtra("Start"));
            str_time.add(data.getStringExtra("Time"));
            str_distance.add(data.getStringExtra("Distance"));
        }
    }
}