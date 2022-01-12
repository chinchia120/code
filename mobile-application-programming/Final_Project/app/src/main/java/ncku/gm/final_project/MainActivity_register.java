package ncku.gm.final_project;

import androidx.appcompat.app.AppCompatActivity;

import android.content.ContentValues;
import android.content.Context;
import android.content.Intent;
import android.database.Cursor;
import android.database.sqlite.SQLiteDatabase;
import android.os.Bundle;
import android.os.Handler;
import android.os.Message;
import android.view.View;
import android.widget.Button;
import android.widget.EditText;

import com.google.android.material.snackbar.Snackbar;

public class MainActivity_register extends AppCompatActivity implements View.OnClickListener {

    SQLiteDatabase db;

    ClientThread mClientThread;
    private Handler mInputHandler;
    private String host = "140.116.47.94";
    private int port = 7070;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main_register);

        db = openOrCreateDatabase("Test_DB", Context.MODE_PRIVATE,null);
        db.execSQL("CREATE TABLE IF NOT EXISTS table01 (_id INTEGER PRIMARY KEY AUTOINCREMENT,name VARCHAR(32),email VARCHAR(32),password VARCHAR(32),phone VARCHAR(32))");

        ((Button)findViewById(R.id.btn_ok_register)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_cancel_register)).setOnClickListener(this);

        mInputHandler = new Handler() {
            @Override
            public void handleMessage(Message msg) {
                if (msg.what == 0) {
                    String[] tmp = msg.obj.toString().split(",");
                    if(tmp[0].matches("user_data")){
                        Cursor cus = db.rawQuery("SELECT * FROM table01",null);
                        int flag=0;
                        cus.moveToFirst();
                        for(int i=0;i<cus.getCount();i++){
                            if(cus.getString(1).matches(tmp[1]) && cus.getString(2).matches(tmp[2]) && cus.getString(3).matches(tmp[3]) && cus.getString(4).matches(tmp[4])){
                                flag=1;
                                break;
                            }else {
                                cus.moveToNext();
                            }
                        }
                        if (flag==0){
                            ContentValues cv = new ContentValues(4);
                            cv.put("name",tmp[1]);
                            cv.put("email",tmp[2]);
                            cv.put("password",tmp[3]);
                            cv.put("phone",tmp[4]);
                            db.insert("table01",null,cv);
                        }
                    }
                }
            }
        };

        mClientThread = new ClientThread(mInputHandler, host, port);
        new Thread(mClientThread).start();
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

                Cursor cus = db.rawQuery("SELECT * FROM table01",null);
                cus.moveToFirst();
                for (int i=0;i<cus.getCount();i++){
                    Message msg = new Message();
                    msg.what = 1;
                    msg.obj = "user_data"+","+
                            cus.getString(1)+","+
                            cus.getString(2)+","+
                            cus.getString(3)+","+
                            cus.getString(4);
                    mClientThread.mOutputHandler.sendMessage(msg);
                    cus.moveToNext();
                }
                finish();
            }
        }else if(view.getId()==R.id.btn_cancel_register){
            finish();
        }
    }
}