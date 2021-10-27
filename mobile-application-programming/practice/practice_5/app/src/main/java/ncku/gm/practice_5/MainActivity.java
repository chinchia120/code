package ncku.gm.practice_5;

import androidx.appcompat.app.AppCompatActivity;

import android.os.Bundle;
import android.widget.Adapter;
import android.widget.ArrayAdapter;
import android.widget.RadioButton;
import android.widget.RadioGroup;
import android.widget.Spinner;
import android.widget.TextView;

import java.util.ArrayList;

public class MainActivity extends AppCompatActivity {

    TextView txv;
    RadioGroup rg_dir,rg_tic;
    Spinner spn_car,spn_seat,spn_go;
    RadioButton rbtn_nor,rbtn_sou;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        txv = findViewById(R.id.txv_msg);

        rg_dir = findViewById(R.id.rg_direction);
        rg_tic = findViewById(R.id.rg_ticket);

        spn_car = findViewById(R.id.spn_car);
        spn_seat = findViewById(R.id.spn_seat);
        spn_go = findViewById(R.id.spn_go_station);

        rbtn_nor = findViewById(R.id.rbtn_north);
        rbtn_sou = findViewById(R.id.rbtn_south);

        String[] car = {"標準車廂","商務車廂"};
        ArrayAdapter<String> ad_car = new ArrayAdapter<>(this, android.R.layout.simple_spinner_dropdown_item,car);
        ((Spinner)findViewById(R.id.spn_car)).setAdapter(ad_car);

        String[] seat = {"靠窗","靠走道"};
        ArrayAdapter<String> ad_seat = new ArrayAdapter<>(this, android.R.layout.simple_spinner_dropdown_item,seat);
        ((Spinner)findViewById(R.id.spn_seat)).setAdapter(ad_seat);

        ArrayList<String> go = new ArrayList<>();
        go.add("台北");
        go.add("桃園");
        go.add("新竹");
        go.add("台中");
        go.add("嘉義");
        go.add("台南");
        go.add("高雄");
        
        ArrayAdapter<String> ad_go = new ArrayAdapter<>(this, android.R.layout.simple_spinner_dropdown_item,go);
        ((Spinner)findViewById(R.id.spn_go_station)).setAdapter(ad_go);

        rg_dir.setOnCheckedChangeListener(new RadioGroup.OnCheckedChangeListener() {
            @Override
            public void onCheckedChanged(RadioGroup radioGroup, int i) {
                if(rbtn_nor.isChecked()){
                    if(go.contains("高雄")){
                        go.remove("台北");
                    }else{
                        go.remove("台北");
                        go.add("高雄");
                    }
                    ((Spinner)findViewById(R.id.spn_go_station)).setAdapter(ad_go);
                }else{
                    if(go.contains("台北")){
                        go.remove("高雄");
                    }else{
                        go.remove("高雄");
                        go.add(0,"台北");
                    }
                    ((Spinner)findViewById(R.id.spn_go_station)).setAdapter(ad_go);
                }
            }
        });

        
    }
}