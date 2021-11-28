package ncku.gm.practice_8;

import androidx.appcompat.app.AppCompatActivity;

import android.content.Intent;
import android.os.Bundle;
import android.view.View;
import android.widget.AdapterView;
import android.widget.ArrayAdapter;
import android.widget.LinearLayout;
import android.widget.ListView;

public class MainActivity4 extends AppCompatActivity implements AdapterView.OnItemClickListener {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main4);

        ((ListView)findViewById(R.id.lsv)).setOnItemClickListener(this);
        String [] str = {"關卡1","關卡2","關卡3"};
        ArrayAdapter<String> ad = new ArrayAdapter<>(this, android.R.layout.simple_list_item_1,str);
        ((ListView)findViewById(R.id.lsv)).setAdapter(ad);
    }

    @Override
    public void onItemClick(AdapterView<?> adapterView, View view, int i, long l) {
        if(i==0){
            startActivity(new Intent(this,MainActivity.class));
        }else if(i==1){
            startActivity(new Intent(this,MainActivity2.class));
        }else if(i==2){
            startActivity(new Intent(this,MainActivity3.class));
        }
    }
}