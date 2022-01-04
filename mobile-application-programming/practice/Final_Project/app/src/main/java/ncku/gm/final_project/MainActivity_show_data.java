package ncku.gm.final_project;

import androidx.appcompat.app.AppCompatActivity;

import android.content.Context;
import android.content.Intent;
import android.database.Cursor;
import android.database.sqlite.SQLiteDatabase;
import android.os.Bundle;
import android.view.View;
import android.widget.AdapterView;
import android.widget.ArrayAdapter;
import android.widget.ImageView;
import android.widget.ListView;
import android.widget.Spinner;
import android.widget.TextView;

import java.io.Serializable;
import java.util.ArrayList;

public class MainActivity_show_data extends AppCompatActivity implements AdapterView.OnItemClickListener , View.OnClickListener {

    ArrayList<String> str_end = new ArrayList<>();
    SQLiteDatabase db;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main_show_data);

        db = openOrCreateDatabase("Test_DB", Context.MODE_PRIVATE,null);
        db.execSQL("CREATE TABLE IF NOT EXISTS table_location (_id INTEGER PRIMARY KEY AUTOINCREMENT,name VARCHAR(32),end_place VARCHAR(32),start VARCHAR(32),time VARCHAR(32))");

        Cursor cus = db.rawQuery("SELECT * FROM table_location",null);
        if(cus.moveToFirst()){
            do{
                str_end.add(cus.getString(2));
            }while (cus.moveToNext());
        }
        ((ListView)findViewById(R.id.lsv_data)).setAdapter(new ArrayAdapter<>(this, android.R.layout.simple_list_item_1,str_end));

        ((ImageView)findViewById(R.id.imv_back_end)).setOnClickListener(this);
        ((ListView)findViewById(R.id.lsv_data)).setOnItemClickListener(this);
    }



    @Override
    public void onItemClick(AdapterView<?> adapterView, View view, int i, long l) {
        Intent it = new Intent(this,MainActivity_show_data_detail.class);
        Bundle bdl = new Bundle();
        bdl.putInt("Cnt",i);
        it.putExtra("Bundle",bdl);
        startActivity(it);
    }

    @Override
    public void onClick(View view) {
        if(view.getId()==R.id.imv_back_end){
            finish();
        }
    }
}