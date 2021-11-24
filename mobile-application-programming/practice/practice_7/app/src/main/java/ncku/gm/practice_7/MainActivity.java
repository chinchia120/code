package ncku.gm.practice_7;

import androidx.annotation.Nullable;
import androidx.appcompat.app.AppCompatActivity;

import android.app.SearchManager;
import android.content.Intent;
import android.net.Uri;
import android.os.Bundle;
import android.view.View;
import android.widget.AdapterView;
import android.widget.ArrayAdapter;
import android.widget.Button;
import android.widget.ListView;

import com.google.android.material.snackbar.Snackbar;

import java.io.Serializable;
import java.util.ArrayList;

public class MainActivity extends AppCompatActivity
        implements AdapterView.OnItemClickListener , AdapterView.OnItemLongClickListener ,
        View.OnClickListener {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        ((ListView)findViewById(R.id.lsv_restaurant)).setOnItemClickListener(this);
        ((ListView)findViewById(R.id.lsv_restaurant)).setOnItemLongClickListener(this);

        ((Button)findViewById(R.id.btn_reset)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_book)).setOnClickListener(this);

        str_res.add("小妞炒飯");
        str_res.add("牛伯麵店");
        str_res.add("余家蔬菜麵");
        str_res.add("麥當勞");
        str_res.add("肯德基");
        str_res.add("新增店家");
        ad_res = new ArrayAdapter<>(this, android.R.layout.simple_list_item_1,str_res);
        ((ListView)findViewById(R.id.lsv_restaurant)).setAdapter(ad_res);
    }

    ArrayList<String> str_res = new ArrayList<>();
    ArrayAdapter<String> ad_res;

    @Override
    public void onItemClick(AdapterView<?> adapterView, View view, int i, long l) {
        if(i==str_res.size()-1){
            startActivityForResult(new Intent(this,MainActivity2.class),000);
        }else{
            startActivity(new Intent(Intent.ACTION_WEB_SEARCH).putExtra(SearchManager.QUERY,((ListView)findViewById(R.id.lsv_restaurant)).getItemAtPosition(i).toString()));
        }
    }

    @Override
    public boolean onItemLongClick(AdapterView<?> adapterView, View view, int i, long l) {
        if(i==str_res.size()-1){
            Snackbar.make(findViewById(R.id.root),"無法刪除",Snackbar.LENGTH_LONG).show();
        }else{
            Snackbar.make(findViewById(R.id.root),"刪除 "+((ListView)findViewById(R.id.lsv_restaurant)).getItemAtPosition(i).toString(),Snackbar.LENGTH_LONG).show();
            str_res.remove(i);
            ad_res.notifyDataSetChanged();
        }
        return true;
    }

    @Override
    public void onClick(View view) {
        ArrayList<String> str = new ArrayList<>();
        for(int i=0;i<str_res.size()-1;i++){
            str.add(i,str_res.get(i));
        }
        if(view.getId()==R.id.btn_book){
            Intent it = new Intent(this,MainActivity3.class);
            Bundle bdl = new Bundle();
            bdl.putSerializable("ArrayList",(Serializable) str);
            it.putExtra("餐廳",bdl);
            startActivityForResult(it,001);
        }else if(view.getId()==R.id.btn_reset){
            str_res.clear();
            str_res.add("小妞炒飯");
            str_res.add("牛伯麵店");
            str_res.add("余家蔬菜麵");
            str_res.add("麥當勞");
            str_res.add("肯德基");
            str_res.add("新增店家");
            ad_res.notifyDataSetChanged();
        }
    }

    @Override
    protected void onActivityResult(int requestCode, int resultCode, @Nullable Intent data) {
        super.onActivityResult(requestCode, resultCode, data);

        if(requestCode==000 && resultCode==-1){
            str_res.add(str_res.size()-1,data.getStringExtra("店名"));
            ad_res.notifyDataSetChanged();
            Snackbar.make(findViewById(R.id.root),"已新增店家",Snackbar.LENGTH_LONG).show();
        }else if(requestCode==001 && resultCode==-1){
            if(data.getIntExtra("確認",0)==-1){
                Snackbar.make(findViewById(R.id.root),"確認訂購",Snackbar.LENGTH_LONG).show();
            }else{
                Snackbar.make(findViewById(R.id.root),"取消訂單",Snackbar.LENGTH_LONG).show();
            }
        }
    }
}