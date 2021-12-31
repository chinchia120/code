package ncku.gm.final_project;

import androidx.appcompat.app.AppCompatActivity;

import android.content.ContentValues;
import android.content.Context;
import android.content.Intent;
import android.database.sqlite.SQLiteDatabase;
import android.os.Bundle;
import android.view.View;
import android.widget.Button;
import android.widget.EditText;

import com.google.android.material.snackbar.Snackbar;

public class MainActivity_register extends AppCompatActivity implements View.OnClickListener {

    SQLiteDatabase db;


    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main_register);

        db = openOrCreateDatabase("Test_DB", Context.MODE_PRIVATE,null);
        db.execSQL("CREATE TABLE IF NOT EXISTS table01 (_id INTEGER PRIMARY KEY AUTOINCREMENT,name VARCHAR(32),email VARCHAR(32),password VARCHAR(32),phone VARCHAR(32))");

        ((Button)findViewById(R.id.btn_ok_register)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_cancel_register)).setOnClickListener(this);
    }

    @Override
    public void onClick(View view) {
        if(view.getId()==R.id.btn_ok_register){
            if(((EditText)findViewById(R.id.edt_enter_name)).getText().toString().matches("") || ((EditText)findViewById(R.id.edt_enter_email_register)).getText().toString().matches("") || ((EditText)findViewById(R.id.edt_enter_password_register)).getText().toString().matches("") || ((EditText)findViewById(R.id.edt_enter_phone)).getText().toString().matches("")){
                Snackbar.make(findViewById(R.id.root_register),"資訊不完整",Snackbar.LENGTH_LONG).show();
            }else{
                ContentValues cv = new ContentValues(4);
                cv.put("name",((EditText)findViewById(R.id.edt_enter_name)).getText().toString());
                cv.put("email",((EditText)findViewById(R.id.edt_enter_email_register)).getText().toString());
                cv.put("password",((EditText)findViewById(R.id.edt_enter_password_register)).getText().toString());
                cv.put("phone",((EditText)findViewById(R.id.edt_enter_phone)).getText().toString());
                db.insert("table01",null,cv);
                finish();
            }
        }else if(view.getId()==R.id.btn_cancel_register){
            finish();
        }
    }
}