package ncku.gm.final_project;

import android.content.ContentValues;
import android.content.Context;
import android.content.DialogInterface;
import android.database.Cursor;
import android.database.sqlite.SQLiteDatabase;
import android.location.Address;
import android.location.Geocoder;
import android.os.Bundle;
import android.os.Handler;
import android.os.Message;
import android.view.View;
import android.view.Menu;

import com.google.android.material.snackbar.Snackbar;
import com.google.android.material.navigation.NavigationView;

import androidx.appcompat.app.AlertDialog;
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
import java.util.Calendar;
import java.util.List;
import java.util.Locale;

public class MainActivity_home extends AppCompatActivity implements OnMapReadyCallback , LocationListener , View.OnClickListener  , DialogInterface.OnClickListener{

    private AppBarConfiguration mAppBarConfiguration;
    private ActivityMainHomeBinding binding;

    private GoogleMap mMap;
    SQLiteDatabase db;

    ClientThread mClientThread;
    private Handler mInputHandler;
    private String host = "140.116.47.94";
    private int port = 7070;

    UserInformation userInformation = new UserInformation();
    LocationData locationData = new LocationData();
    /*
    //對話框debug
    int diaread = 0;
     */


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
        lm.requestLocationUpdates(LocationManager.NETWORK_PROVIDER, 1000, 5, this);
        lm.requestLocationUpdates(LocationManager.GPS_PROVIDER, 1000, 5, this);

        SupportMapFragment smf = (SupportMapFragment) getSupportFragmentManager().findFragmentById(R.id.map);
        smf.getMapAsync(this);

        ((Button)findViewById(R.id.btn_new_data)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_show_data)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_user_data)).setOnClickListener(this);

        db = openOrCreateDatabase("Test_DB", Context.MODE_PRIVATE,null);
        db.execSQL("CREATE TABLE IF NOT EXISTS table_location (_id INTEGER PRIMARY KEY AUTOINCREMENT,name VARCHAR(32),end_place VARCHAR(32),start VARCHAR(32),time VARCHAR(32))");

        mInputHandler = new Handler()
        {
            @Override
            public void handleMessage(Message msg)
            {
                if (msg.what == 0)
                {
                    String[] tmp = msg.obj.toString().split(",");
                    //收到新訊息
                    if(tmp[0].matches("message")){
                        if (!tmp[1].matches(userInformation.getUser_name())){
                            Snackbar.make(findViewById(R.id.drawer_layout),"Get Message",Snackbar.LENGTH_LONG).show();
                        }
                    }else if(tmp[0].matches("together")){
                        if(tmp[1].matches(userInformation.getUser_name())){

                            locationData.setLat_start(Double.valueOf(tmp[2]));
                            locationData.setLon_start(Double.valueOf(tmp[3]));
                            locationData.setLat_end(Double.valueOf(tmp[4]));
                            locationData.setLon_end(Double.valueOf(tmp[5]));
                            check_together();
                        }
                    }else if(tmp[0].matches("delete")){
                        Cursor cus = db.rawQuery("SELECT * FROM table_location",null);
                        cus.moveToFirst();

                        for (int i=0;i<cus.getCount();i++){
                            if(tmp[1].matches(cus.getString(1)) && tmp[2].matches(cus.getString(2)) && tmp[3].matches(cus.getString(3)) && tmp[4].matches(cus.getString(4))){
                                db.delete("table_location","_id="+cus.getString(0),null);
                                break;
                            }
                            cus.moveToNext();
                        }
                    }
                }
            }
        };

        mClientThread = new ClientThread(mInputHandler, host, port);
        new Thread(mClientThread).start();
    }

        private void check_together(){
            AlertDialog.Builder b = new AlertDialog.Builder(this);
            b.setIcon(android.R.drawable.presence_away);
            b.setMessage("有人要與您共乘");
            b.setPositiveButton("查看",this);
            b.setCancelable(false);
            b.show();
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
        mMap.clear();
    }

    @Override
    public void onLocationChanged(@NonNull Location location) {
        if (location != null){
            if (mMap != null){
                mMap.animateCamera(CameraUpdateFactory.newLatLngZoom(new LatLng(location.getLatitude(),location.getLongitude()),7f));
                mMap.clear();
                mMap.addMarker(new MarkerOptions().position(new LatLng(location.getLatitude(),location.getLongitude())).title("目前位置"));
                
                Cursor cus = db.rawQuery("SELECT * FROM table_location",null);
                if(cus.moveToFirst()){
                    do{
                        Geocoder geo = new Geocoder(this, Locale.TRADITIONAL_CHINESE);
                        try {
                            List<Address> list_start = geo.getFromLocationName(cus.getString(3),1);
                            mMap.addMarker(new MarkerOptions().position(new LatLng(list_start.get(0).getLatitude(),list_start.get(0).getLongitude())).title("出發位置"));
                        }catch (Exception e){
                            e.printStackTrace();
                        }
                    }while (cus.moveToNext());
                }
            }
        }

        Cursor cus = db.rawQuery("SELECT * FROM table_location",null);
        cus.moveToFirst();
        for(int i=0;i<cus.getCount();i++){
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
            if (nowmon == startmonth && nowdate == startdate && lefttime == 0){
                db.delete("table_location","_id="+cus.getString(0),null);
            }
            cus.moveToNext();
        }
    }

    @Override
    public void onClick(View view) {
        if(view.getId()==R.id.btn_new_data){
            startActivity(new Intent(this,MainActivity_new_data.class));
        }else if(view.getId()==R.id.btn_show_data){
            startActivity(new Intent(this,MainActivity_show_data.class));
        }else if(view.getId()==R.id.btn_user_data){
            startActivity(new Intent(this,MainActivity_user_data.class));
        }
    }

    @Override
    public void onClick(DialogInterface dialogInterface, int i) {
        if(i==-1){
            startActivity(new Intent(this,MainActivity_together.class));
        }
    }
}