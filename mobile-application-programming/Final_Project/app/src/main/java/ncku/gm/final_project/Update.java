package ncku.gm.final_project;

import android.content.DialogInterface;
import android.os.Handler;
import androidx.annotation.NonNull;
import androidx.appcompat.app.AlertDialog;
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

import com.google.android.material.snackbar.Snackbar;

import java.util.HashMap;
import java.util.Random;
import java.util.jar.Manifest;

public class Update {
    ClientThread mClientThread;
    private Handler mInputHandler;
    private String host = "140.116.47.94";
    private int port = 7070;

    public void check(){
        mInputHandler = new Handler()
        {
            @Override
            public void handleMessage(Message msg)
            {
                if (msg.what == 0)
                {
                    String[] tmp = msg.obj.toString().split(",");
                    //私訊
                    if(tmp[0].matches("together")){
                        AlertDialog.Builder bdr = new AlertDialog.Builder(null);
                        bdr.setTitle("Together");
                        bdr.show();
                    }
                }
            }
        };

        mClientThread = new ClientThread(mInputHandler, host, port);
        new Thread(mClientThread).start();
    }


}
