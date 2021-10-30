package ncku.gm.practice_5;

import androidx.appcompat.app.AppCompatActivity;

import android.os.Bundle;
import android.view.View;
import android.widget.AdapterView;
import android.widget.ArrayAdapter;
import android.widget.ListView;
import android.widget.RadioButton;
import android.widget.RadioGroup;
import android.widget.Spinner;
import android.widget.TextView;

import java.util.ArrayList;

public class MainActivity extends AppCompatActivity
        implements RadioGroup.OnCheckedChangeListener, AdapterView.OnItemSelectedListener, AdapterView.OnItemClickListener{

    TextView txv;
    RadioGroup rg_dir,rg_tic;
    Spinner spn_car,spn_seat,spn_go;
    RadioButton rbtn_nor,rbtn_sou,rbtn_adu,rntn_chi;
    ListView lsv_stop;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        txv = findViewById(R.id.txv_msg);

        rg_dir = findViewById(R.id.rg_direction);
        rg_dir.setOnCheckedChangeListener(this);
        rg_tic = findViewById(R.id.rg_ticket);

        spn_car = findViewById(R.id.spn_car);
        spn_seat = findViewById(R.id.spn_seat);
        spn_go = findViewById(R.id.spn_go_station);
        spn_go.setOnItemSelectedListener(this);

        rbtn_nor = findViewById(R.id.rbtn_north);
        rbtn_sou = findViewById(R.id.rbtn_south);
        rbtn_adu = findViewById(R.id.rbtn_adult);
        rntn_chi = findViewById(R.id.rbtn_child);

        lsv_stop = findViewById(R.id.lsv_stop_station);
        lsv_stop.setOnItemClickListener(this);

        String[] car = {"標準車廂","自由座"};
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

        ArrayList<String> stop = new ArrayList<>();
        stop.add("台北");
        stop.add("桃園");
        stop.add("新竹");
        stop.add("台中");
        stop.add("嘉義");
        stop.add("台南");
        stop.add("高雄");
        ArrayAdapter<String> ad_stop = new ArrayAdapter<>(this, android.R.layout.simple_list_item_1,stop);
        ((ListView)findViewById(R.id.lsv_stop_station)).setAdapter(ad_stop);
    }

    @Override
    public void onCheckedChanged(RadioGroup radioGroup, int i) {
        ArrayList<String> go = new ArrayList<>();
        if(rbtn_nor.isChecked()){
            go.add("桃園");
            go.add("新竹");
            go.add("台中");
            go.add("嘉義");
            go.add("台南");
            go.add("高雄");
        }else{
            go.add("台北");
            go.add("桃園");
            go.add("新竹");
            go.add("台中");
            go.add("嘉義");
            go.add("台南");
        }
        ArrayAdapter<String> ad_go = new ArrayAdapter<>(this, android.R.layout.simple_spinner_dropdown_item,go);
        ((Spinner)findViewById(R.id.spn_go_station)).setAdapter(ad_go);
    }

    @Override
    public void onItemSelected(AdapterView<?> adapterView, View view, int position, long id) {
        ArrayList<String> stop = new ArrayList<>();
        stop.add("台北");
        stop.add("桃園");
        stop.add("新竹");
        stop.add("台中");
        stop.add("嘉義");
        stop.add("台南");
        stop.add("高雄");
        if(rbtn_nor.isChecked()){
            for(int i=position+1;i<7;i++){
                stop.remove(position+1);
            }
        }else{
            for(int i=0;i<position+1;i++){
                stop.remove(0);
            }
        }
        ArrayAdapter<String> ad_stop = new ArrayAdapter<>(this, android.R.layout.simple_list_item_1,stop);
        ((ListView)findViewById(R.id.lsv_stop_station)).setAdapter(ad_stop);
    }

    @Override
    public void onNothingSelected(AdapterView<?> adapterView) {

    }

    @Override
    public void onItemClick(AdapterView<?> adapterView, View view, int position, long id) {
        String str_show = "";
        String str_go = spn_go.getSelectedItem().toString();
        int idx_go = spn_go.getSelectedItemPosition();
        int[] ticket_nor = {160,130,410,380,280,140};
        int[] ticket_free = {155,125,395,365,270,135};
        int tmp_go = 0 , tmp_stop = 0;
        String[] station = {"台北","桃園","新竹","台中","嘉義","台南","高雄"};
        if(rbtn_nor.isChecked()){
            str_show += "北上" + str_go + "到" + station[position] +"\n";
            for(int i=0;i<7;i++){
                if(str_go==station[i]){
                    tmp_go=i;
                    tmp_stop=position;
                }
            }
        }else{
            str_show += "南下" + str_go + "到" + station[position+idx_go+1]+"\n";
            for(int i=0;i<7;i++){
                if(str_go==station[i]){
                    tmp_go=i;
                    tmp_stop=position+idx_go+1;
                }
            }
        }
        String str_car = spn_car.getSelectedItem().toString();
        int cnt = 0;
        int idx_car =spn_car.getSelectedItemPosition();
        if(idx_car==0){
            for(int i=Math.min(tmp_go,tmp_stop);i<Math.max(tmp_go,tmp_stop);i++){
                cnt += ticket_nor[i];
            }
            if(Math.max(tmp_go,tmp_stop)>4){
                cnt -= 10;
            }
            if((tmp_go==4 && tmp_stop==5) || (tmp_go==5 && tmp_stop==4)){
                cnt=280;
            }else if((tmp_go==5 && tmp_stop==6) || (tmp_go==6 && tmp_stop==5)){
                cnt=140;
            }
        }else{
            for(int i=Math.min(tmp_go,tmp_stop);i<Math.max(tmp_go,tmp_stop);i++){
                cnt += ticket_free[i];
            }
        }
        String str_seat = spn_seat.getSelectedItem().toString();
        str_show += str_car + "的" + str_seat + "車位\n";
        if(rbtn_adu.isChecked()){
            str_show += "全票票價一共" + cnt + "元";
        }else{
            str_show += "優惠票票價一共" + String.format("%.0f",cnt*0.5) + "元" ;
        }
        txv.setText(str_show+"");
    }
}