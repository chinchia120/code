package ncku.gm.final_project;

import androidx.annotation.Nullable;
import androidx.appcompat.app.AppCompatActivity;

import android.content.ContentValues;
import android.content.Context;
import android.content.Intent;
import android.database.Cursor;
import android.database.sqlite.SQLiteDatabase;
import android.os.Bundle;
import android.view.View;
import android.widget.Button;
import android.widget.EditText;
import android.widget.ListView;
import android.widget.TextView;

import com.google.android.material.snackbar.Snackbar;

public class MainActivity extends AppCompatActivity implements View.OnClickListener {

    SQLiteDatabase db;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        startActivity(new Intent(this,MainActivity_home.class));

        db = openOrCreateDatabase("Test_DB", Context.MODE_PRIVATE,null);
        db.execSQL("CREATE TABLE IF NOT EXISTS table01 (_id INTEGER PRIMARY KEY AUTOINCREMENT,name VARCHAR(32),email VARCHAR(32),password VARCHAR(32),phone VARCHAR(32))");
        //query();

        ((Button)findViewById(R.id.btn_login)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_register)).setOnClickListener(this);
    }

    @Override
    public void onClick(View view) {
        if(view.getId()==R.id.btn_login){
            Cursor cus = db.rawQuery("SELECT * FROM table01",null);
            if(cus.moveToFirst()){
                do {
                    if(cus.getString(2).equals(((EditText)findViewById(R.id.edt_enter_email)).getText().toString()) && cus.getString(3).equals(((EditText)findViewById(R.id.edt_enter_password)).getText().toString())){
                        startActivity(new Intent(this,MainActivity_home.class));
                    }
                }while (cus.moveToNext());
            }
            Snackbar.make(findViewById(R.id.root_login),"登入失敗",Snackbar.LENGTH_LONG).show();
        }else if(view.getId()==R.id.btn_register){
            startActivityForResult(new Intent(this,MainActivity_register.class),001);
        }
    }

    @Override
    protected void onActivityResult(int requestCode, int resultCode, @Nullable Intent data) {
        super.onActivityResult(requestCode, resultCode, data);

        if(requestCode==001 && resultCode==RESULT_OK){
            Snackbar.make(findViewById(R.id.root_login),"註冊成功",Snackbar.LENGTH_LONG).show();
            ContentValues cv = new ContentValues(4);
            cv.put("name",data.getStringExtra("Name"));
            cv.put("email",data.getStringExtra("Email"));
            cv.put("password",data.getStringExtra("Password"));
            cv.put("phone",data.getStringExtra("Phone"));
            db.insert("table01",null,cv);
            //query();
        }
    }

    public void query(){
        ((TextView)findViewById(R.id.textView)).setText(db.getPath()+"\n"+db.getPageSize()+" Bytes\n"+db.getMaximumSize()+" Bytes\n");
        Cursor cus = db.rawQuery("SELECT * FROM table01",null);
        ((TextView)findViewById(R.id.textView)).append(cus.getCount()+"\n");
        if(cus.moveToFirst()){
            do{
                ((TextView)findViewById(R.id.textView)).append(cus.getString(1)+"\t"+cus.getString(2)+"\t"+cus.getString(3)+"\t"+cus.getString(4)+"\n");
            }while (cus.moveToNext());
        }
    }
}