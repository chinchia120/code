package ncku.gm.final_project;

import androidx.appcompat.app.AppCompatActivity;

import android.os.Bundle;
import android.view.View;
import android.widget.Button;
import android.widget.TextView;

import java.util.ArrayList;

public class MainActivity_show_data_detail extends AppCompatActivity implements View.OnClickListener {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main_show_data_detail);

        ((TextView)findViewById(R.id.txv_show_end)).setText((String) getIntent().getBundleExtra("Bundle").getSerializable("End"));
        ((TextView)findViewById(R.id.txv_show_start)).setText("出發地 : "+(String) getIntent().getBundleExtra("Bundle").getSerializable("Start"));
        ((TextView)findViewById(R.id.txv_show_distance)).setText("距離 : "+(String) getIntent().getBundleExtra("Bundle").getSerializable("Distance"));
        ((TextView)findViewById(R.id.txv_show_time)).setText("時間 : "+(String) getIntent().getBundleExtra("Bundle").getSerializable("Time"));


        ((Button)findViewById(R.id.btn_back_data_detail)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_chat)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_sure)).setOnClickListener(this);
    }

    @Override
    public void onClick(View view) {
        if (view.getId()==R.id.btn_back_data_detail){
            finish();
        }else if(view.getId()==R.id.btn_chat){

        }else if(view.getId()==R.id.btn_sure){

        }
    }
}