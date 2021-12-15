package ncku.gm.p1215;

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
public class MainActivity extends AppCompatActivity implements LocationListener{

    private EditText etMain;
    private Button btnMain;
    private TextView tvMain;
    ClientThread mClientThread;
    private Handler mInputHandler;
    private String host = "140.116.47.94";
    private int port = 7070;
    private HashMap<String, String> userLocation;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        Random r = new Random();
        etMain = findViewById(R.id.et_main);
        etMain.setText(r.nextInt(10000)+"");
        btnMain = findViewById(R.id.btn_main);
        tvMain = findViewById(R.id.tv_main);
        userLocation = new HashMap<String, String>();

        LocationManager lm = (LocationManager) getSystemService(LOCATION_SERVICE);
        if (ActivityCompat.checkSelfPermission(this, android.Manifest.permission.ACCESS_FINE_LOCATION) != PackageManager.PERMISSION_GRANTED && ActivityCompat.checkSelfPermission(this, android.Manifest.permission.ACCESS_COARSE_LOCATION) != PackageManager.PERMISSION_GRANTED) {
            // TODO: Consider calling
            //    ActivityCompat#requestPermissions
            // here to request the missing permissions, and then overriding
            //   public void onRequestPermissionsResult(int requestCode, String[] permissions,
            //                                          int[] grantResults)
            // to handle the case where the user grants the permission. See the documentation
            // for ActivityCompat#requestPermissions for more details.
            ActivityCompat.requestPermissions(this, new String[]{android.Manifest.permission.ACCESS_FINE_LOCATION, android.Manifest.permission.ACCESS_COARSE_LOCATION}, 100);
            return;
        }
        lm.requestLocationUpdates(LocationManager.GPS_PROVIDER, 1000, 5, this);
        lm.requestLocationUpdates(LocationManager.NETWORK_PROVIDER, 1000, 5, this);
        mInputHandler = new Handler()
        {
            @Override
            public void handleMessage(Message msg)
            {
                if (msg.what == 0)
                {
                    String[] tmp = msg.obj.toString().split(",");
                    userLocation.put(tmp[0], tmp[1]+","+tmp[2]);
                    tvMain.setText(userLocation.toString());
                }
            }
        };

        btnMain.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                try {
                    //Toast.makeText(this, etMain.getText().toString()+","+location.getLatitude()+","+location.getLongitude(), Toast.LENGTH_SHORT).show();
                    Message msg = new Message();
                    msg.what = 1;
                    msg.obj = etMain.getText().toString()+",0,0";
                    mClientThread.mOutputHandler.sendMessage(msg);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        });

        mClientThread = new ClientThread(mInputHandler, host, port);
        new Thread(mClientThread).start();
    }

    @Override
    public void onLocationChanged(@NonNull Location location) {
        try {
            //Toast.makeText(this, etMain.getText().toString()+","+location.getLatitude()+","+location.getLongitude(), Toast.LENGTH_SHORT).show();
            Message msg = new Message();
            msg.what = 1;
            msg.obj = etMain.getText().toString()+","+location.getLatitude()+","+location.getLongitude();
            mClientThread.mOutputHandler.sendMessage(msg);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}