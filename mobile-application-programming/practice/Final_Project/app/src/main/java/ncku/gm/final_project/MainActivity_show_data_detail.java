package ncku.gm.final_project;

import androidx.appcompat.app.AppCompatActivity;

import android.content.Context;
import android.content.Intent;
import android.database.Cursor;
import android.database.sqlite.SQLiteDatabase;
import android.os.Bundle;
import android.view.View;
import android.widget.Button;
import android.widget.EditText;
import android.widget.TextView;

import com.google.android.material.snackbar.Snackbar;

import java.util.ArrayList;

public class MainActivity_show_data_detail extends AppCompatActivity implements View.OnClickListener {

    SQLiteDatabase db;
    UserInformation userInformation = new UserInformation();
    Cursor cus;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main_show_data_detail);

        ((Button)findViewById(R.id.btn_back_data_detail)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_chat)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_sure)).setOnClickListener(this);

        db = openOrCreateDatabase("Test_DB", Context.MODE_PRIVATE,null);
        db.execSQL("CREATE TABLE IF NOT EXISTS table02 (_id INTEGER PRIMARY KEY AUTOINCREMENT,name VARCHAR(32),end_place VARCHAR(32),start VARCHAR(32),distance VARCHAR(32),time VARCHAR(32))");
        cus = db.rawQuery("SELECT * FROM table02",null);
        if(cus.moveToFirst()){
            int i = getIntent().getBundleExtra("Bundle").getInt("Cnt");
            while (i>0){
                cus.moveToNext();
                i--;
            }
        }

        ((TextView)findViewById(R.id.txv_show_name)).setText("發起人 : "+cus.getString(1));
        ((TextView)findViewById(R.id.txv_show_end)).setText(cus.getString(2));
        ((TextView)findViewById(R.id.txv_show_start)).setText("出發處 : "+cus.getString(3));
        ((TextView)findViewById(R.id.txv_show_distance)).setText("距離 : "+cus.getString(4));
        ((TextView)findViewById(R.id.txv_show_time)).setText("時間 : "+cus.getString(5));
    }

    @Override
    public void onClick(View view) {
        if (view.getId()==R.id.btn_back_data_detail){
            finish();
        }else if(view.getId()==R.id.btn_chat){

        }else if(view.getId()==R.id.btn_sure){
            if(cus.getString(1).matches(userInformation.getUser_name())){
                Snackbar.make(findViewById(R.id.root_show_data_detail),"這是自己發起的共乘喔",Snackbar.LENGTH_LONG).show();
            }else {

            }
        }
    }
}