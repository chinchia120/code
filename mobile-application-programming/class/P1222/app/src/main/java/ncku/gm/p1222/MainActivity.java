package ncku.gm.p1222;

import androidx.appcompat.app.AppCompatActivity;

import android.content.ContentValues;
import android.content.Context;
import android.database.Cursor;
import android.database.sqlite.SQLiteDatabase;
import android.os.Bundle;
import android.view.View;
import android.widget.AdapterView;
import android.widget.EditText;
import android.widget.ListView;
import android.widget.SimpleCursorAdapter;
import android.widget.TextView;

public class MainActivity extends AppCompatActivity implements AdapterView.OnItemClickListener {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        db = openOrCreateDatabase("Test_DB", Context.MODE_PRIVATE,null);
        //String cmd = "CREATE TABLE IF NOT EXISTS table_1 (name VARCHAR(32),phone VARCHAR(32),email VARCHAR(32))";
        //String cmd = "CREATE TABLE IF NOT EXISTS table_2 (name VARCHAR(32) PRIMARY KEY,phone VARCHAR(32),email VARCHAR(32))";
        String cmd = "CREATE TABLE IF NOT EXISTS table_3 (_id INTEGER PRIMARY KEY AUTOINCREMENT,name VARCHAR(32),phone VARCHAR(32),email VARCHAR(32))";
        db.execSQL(cmd);

        //add_Data("AAA","0987654321","AAA@gmail.com");
        query();

        ((ListView)findViewById(R.id.lsv_show)).setOnItemClickListener(this);
    }

    SQLiteDatabase db;

    private void add_Data(String name , String phone , String email){
        ContentValues cv = new ContentValues(3);
        cv.put("name",name);
        cv.put("phone",phone);
        cv.put("email",email);
        //db.insert("table_1",null,cv);
        //db.insert("table_2",null,cv);
        db.insert("table_3",null,cv);
        query();
    }

    private void query(){
        //((TextView)findViewById(R.id.txv_show)).setText(db.getPath()+"\n"+db.getPageSize()+" Bytes\n"+db.getMaximumSize()+" Bytes\n");
        //Cursor cus = db.rawQuery("SELECT * FROM table_1",null);
        //Cursor cus = db.rawQuery("SELECT * FROM table_2",null);
        Cursor cus = db.rawQuery("SELECT * FROM table_3",null);

        SimpleCursorAdapter sca = new SimpleCursorAdapter(this,R.layout.item,cus,new String[]{"name","phone","email"},new int[]{R.id.txv_name,R.id.txv_phone,R.id.txv_email},0);
        ((ListView)findViewById(R.id.lsv_show)).setAdapter(sca);

        //((TextView)findViewById(R.id.txv_show)).append(cus.getCount()+"\n");
         /*
        if(cus.moveToFirst()){
            do{
                String name = cus.getString(1);
                String phone = cus.getString(2);
                String email = cus.getString(3);
                ((TextView)findViewById(R.id.txv_show)).append(name+"\t"+phone+"\t"+email+"\n");
            }while (cus.moveToNext());
        }
         */
    }

    public void new_Data(View v){
        String name = ((EditText)findViewById(R.id.edt_name)).getText().toString();
        String phone = ((EditText)findViewById(R.id.edt_phone)).getText().toString();
        String email = ((EditText)findViewById(R.id.edt_email)).getText().toString();
        add_Data(name,phone,email);
        query();
    }

    public void set_Data(View v){
        String name = ((EditText)findViewById(R.id.edt_name)).getText().toString();
        String phone = ((EditText)findViewById(R.id.edt_phone)).getText().toString();
        String email = ((EditText)findViewById(R.id.edt_email)).getText().toString();
        update(name,phone,email,_id);
        query();
    }

    private void update(String name , String phone ,String email , int id){
        ContentValues cv = new ContentValues(3);
        cv.put("name",name);
        cv.put("phone",phone);
        cv.put("email",email);
        db.update("table_3",cv,"_id="+id,null);
    }

    public void delete_Data(View v){
        db.delete("table_3","_id="+_id,null);
        ((EditText)findViewById(R.id.edt_name)).setText("");
        ((EditText)findViewById(R.id.edt_phone)).setText("");
        ((EditText)findViewById(R.id.edt_email)).setText("");
        query();
    }

    int _id=0;
    @Override
    public void onItemClick(AdapterView<?> adapterView, View view, int i, long l) {
        Cursor cus = db.rawQuery("SELECT * FROM table_3",null);
        if(cus.moveToFirst()){
            while (i>0){
                cus.moveToNext();
                i--;
            }
        }
        _id = cus.getInt(0);
        String name = cus.getString(1);
        String phone = cus.getString(2);
        String email = cus.getString(3);
        ((EditText)findViewById(R.id.edt_name)).setText(name);
        ((EditText)findViewById(R.id.edt_phone)).setText(phone);
        ((EditText)findViewById(R.id.edt_email)).setText(email);
    }
}