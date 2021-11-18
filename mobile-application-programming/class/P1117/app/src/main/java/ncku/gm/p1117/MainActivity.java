package ncku.gm.p1117;

import androidx.annotation.Nullable;
import androidx.appcompat.app.AppCompatActivity;

import android.content.Intent;
import android.net.Uri;
import android.os.Bundle;
import android.view.View;
import android.widget.AdapterView;
import android.widget.ArrayAdapter;
import android.widget.Button;
import android.widget.ListView;

public class MainActivity extends AppCompatActivity
        implements View.OnClickListener , AdapterView.OnItemClickListener {

    ArrayAdapter<String> ad_note;
    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        ((Button)findViewById(R.id.btn_call_act2)).setOnClickListener(this);

        ((ListView)findViewById(R.id.lsv_note)).setOnItemClickListener(this);

        ad_note = new ArrayAdapter<>(this, android.R.layout.simple_list_item_1,str_note);
        ((ListView)findViewById(R.id.lsv_note)).setAdapter(ad_note);

    }

    String [] str_note = {"1. 按一下可以編輯","2. 長按可以清除","3.","4.","5.","6."};

    @Override
    public void onClick(View view) {
        if(view.getId()==R.id.btn_call_act2){
            Intent it = new Intent();
            it.setAction(Intent.ACTION_VIEW);
            //Uri uri = Uri.parse("http://flag.com.tw");
            //it.setData(uri);
            //new Intent().setAction(Intent.ACTION_VIEW);
            //new Intent().setData((Uri.parse("http://flag.com.tw")));
            it.setData(Uri.parse("http://flag.com.tw"));
            //startActivity(new Intent(this,MainActivity2.class));
            startActivity(it);
        }
    }

    @Override
    public void onItemClick(AdapterView<?> adapterView, View view, int i, long l) {
        if(adapterView.getId()==R.id.lsv_note){
            Intent it = new Intent();
            it.setClass(this,MainActivity2.class);
            it.putExtra("編號",i+1);
            if(str_note[i].length()>3){
                it.putExtra("內容",str_note[i].substring(3));
            }else{
                it.putExtra("內容","");
            }
            startActivityForResult(it,000);
        }
    }

    @Override
    protected void onActivityResult(int requestCode, int resultCode, @Nullable Intent data) {
        super.onActivityResult(requestCode, resultCode, data);
        if(requestCode==000 && resultCode==-1){
            int id = data.getIntExtra("編號",0);
            String str_text = data.getStringExtra("內容");
            String str_show = id + ". " + str_text;

            str_note [id-1] = str_show;
            //ArrayAdapter<String> ad_note = new ArrayAdapter<>(this, android.R.layout.simple_list_item_1,str_note);
            //((ListView)findViewById(R.id.lsv_note)).setAdapter(ad_note);
            ad_note.notifyDataSetChanged();
        }
    }
}