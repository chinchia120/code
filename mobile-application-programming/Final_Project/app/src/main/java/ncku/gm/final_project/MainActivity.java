package ncku.gm.final_project;

import androidx.annotation.Nullable;
import androidx.appcompat.app.AppCompatActivity;

import android.content.ContentValues;
import android.content.Context;
import android.content.Intent;
import android.database.Cursor;
import android.database.sqlite.SQLiteDatabase;
import android.os.Bundle;
import android.service.autofill.UserData;
import android.view.View;
import android.widget.Button;
import android.widget.EditText;
import android.widget.ListView;
import android.widget.TextView;

import com.google.android.gms.maps.model.LatLng;
import com.google.android.gms.maps.model.MarkerOptions;
import com.google.android.material.snackbar.Snackbar;

import java.io.Serializable;

public class MainActivity extends AppCompatActivity implements View.OnClickListener {

    SQLiteDatabase db,db_location;
    UserInformation userInformation = new UserInformation();

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        //startActivity(new Intent(this,MainActivity_home.class));

        db = openOrCreateDatabase("Test_DB", Context.MODE_PRIVATE,null);
        db.execSQL("CREATE TABLE IF NOT EXISTS table01 (_id INTEGER PRIMARY KEY AUTOINCREMENT,name VARCHAR(32),email VARCHAR(32),password VARCHAR(32),phone VARCHAR(32))");
        //query();
        //show_location();

        ((Button)findViewById(R.id.btn_login)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_register)).setOnClickListener(this);
    }

    @Override
    public void onClick(View view) {
        if(view.getId()==R.id.btn_login){
            //query();
            //show_location();
            Cursor cus = db.rawQuery("SELECT * FROM table01",null);
            if(cus.moveToFirst()){
                do {
                    if(cus.getString(2).equals(((EditText)findViewById(R.id.edt_enter_email)).getText().toString()) && cus.getString(3).equals(((EditText)findViewById(R.id.edt_enter_password)).getText().toString())){
                        userInformation.setUser_name(cus.getString(1));
                        userInformation.setUser_email(cus.getString(2));
                        userInformation.setUser_password(cus.getString(3));
                        userInformation.setUser_phone(cus.getString(4));
                        startActivity(new Intent(this,MainActivity_home.class));
                        return;
                    }
                }while (cus.moveToNext());
            }
            Snackbar.make(findViewById(R.id.root_login),"登入失敗",Snackbar.LENGTH_LONG).show();
        }else if(view.getId()==R.id.btn_register){
            startActivity(new Intent(this,MainActivity_register.class));
        }
    }

    private void query(){
        ((TextView)findViewById(R.id.textView)).setText(db.getPath()+"\n"+db.getPageSize()+" Bytes\n"+db.getMaximumSize()+" Bytes\n");
        Cursor cus = db.rawQuery("SELECT * FROM table01",null);
        ((TextView)findViewById(R.id.textView)).append(cus.getCount()+"\n");
/*
        cus.moveToFirst();
        for(int i=0;i<cus.getCount();i++){
            db.delete("table01","_id="+cus.getString(0),null);
            cus.moveToNext();
        }
 */

        if(cus.moveToFirst()){
            do{
                ((TextView)findViewById(R.id.textView)).append(cus.getString(0)+"\t"+cus.getString(1)+"\t"+cus.getString(2)+"\t"+cus.getString(3)+"\t"+cus.getString(4)+"\n");
            }while (cus.moveToNext());
        }
    }

    private void show_location(){
        ((TextView)findViewById(R.id.textView)).setText(db.getPath()+"\n"+db.getPageSize()+" Bytes\n"+db.getMaximumSize()+" Bytes\n");
        db_location = openOrCreateDatabase("Test_DB", Context.MODE_PRIVATE,null);
        db_location.execSQL("CREATE TABLE IF NOT EXISTS table_location (_id INTEGER PRIMARY KEY AUTOINCREMENT,name VARCHAR(32),end_place VARCHAR(32),start VARCHAR(32),time VARCHAR(32))");
        Cursor cus_location = db_location.rawQuery("SELECT * FROM table_location",null);
        ((TextView)findViewById(R.id.textView)).append(cus_location.getCount()+"\n");
        if(cus_location.moveToFirst()){
            do{
                ((TextView)findViewById(R.id.textView)).append(cus_location.getString(1)+"\t"+cus_location.getString(2)+"\t"+cus_location.getString(3)+"\t"+cus_location.getString(4)+"\n");
            }while (cus_location.moveToNext());
        }
    }
}