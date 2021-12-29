package ncku.gm.final_project;

import androidx.appcompat.app.AppCompatActivity;

import android.content.Intent;
import android.os.Bundle;
import android.view.View;
import android.widget.AdapterView;
import android.widget.ArrayAdapter;
import android.widget.ImageView;
import android.widget.ListView;
import android.widget.Spinner;

import java.io.Serializable;
import java.util.ArrayList;

public class MainActivity_show_data extends AppCompatActivity implements AdapterView.OnItemClickListener , View.OnClickListener {

    ArrayList<String> str_end = new ArrayList<>();
    ArrayList<String> str_start = new ArrayList<>();
    ArrayList<String> str_distance = new ArrayList<>();
    ArrayList<String> str_time = new ArrayList<>();

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main_show_data);

        str_end = (ArrayList<String>) getIntent().getBundleExtra("Bundle").getSerializable("End");
        str_start = (ArrayList<String>) getIntent().getBundleExtra("Bundle").getSerializable("Start");
        str_time = (ArrayList<String>) getIntent().getBundleExtra("Bundle").getSerializable("Time");
        str_distance = (ArrayList<String>) getIntent().getBundleExtra("Bundle").getSerializable("Distance");
        ((ListView)findViewById(R.id.lsv_data)).setAdapter(new ArrayAdapter<>(this, android.R.layout.simple_list_item_1,str_end));

        ((ImageView)findViewById(R.id.imv_back_end)).setOnClickListener(this);
        ((ListView)findViewById(R.id.lsv_data)).setOnItemClickListener(this);
    }



    @Override
    public void onItemClick(AdapterView<?> adapterView, View view, int i, long l) {
        Intent it = new Intent(this,MainActivity_show_data_detail.class);
        Bundle bdl = new Bundle();
        bdl.putString("End",str_end.get(i));
        bdl.putString("Start",str_start.get(i));
        bdl.putString("Time",str_time.get(i));
        bdl.putString("Distance",str_distance.get(i));
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