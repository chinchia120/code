package ncku.gm.practice_6;

import androidx.appcompat.app.AppCompatActivity;

import android.app.AlertDialog;
import android.app.DatePickerDialog;
import android.app.TimePickerDialog;
import android.content.DialogInterface;
import android.os.Bundle;
import android.view.View;
import android.widget.AdapterView;
import android.widget.ArrayAdapter;
import android.widget.Button;
import android.widget.CompoundButton;
import android.widget.DatePicker;
import android.widget.ListView;
import android.widget.Spinner;
import android.widget.TextView;
import android.widget.TimePicker;
import android.widget.Toast;

import com.google.android.material.snackbar.Snackbar;

import java.util.ArrayList;

public class MainActivity extends AppCompatActivity
        implements View.OnClickListener , DatePickerDialog.OnDateSetListener ,
        TimePickerDialog.OnTimeSetListener , DialogInterface.OnClickListener ,
        AdapterView.OnItemSelectedListener , AdapterView.OnItemClickListener  ,
        CompoundButton.OnCheckedChangeListener {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        ((Button)findViewById(R.id.btn_book)).setOnClickListener(this);

        ((Spinner)findViewById(R.id.spn_league)).setOnItemSelectedListener(this);
        ((Spinner)findViewById(R.id.spn_area)).setOnItemSelectedListener(this);
        ((Spinner)findViewById(R.id.spn_seat)).setOnItemSelectedListener(this);

        ((ListView)findViewById(R.id.lsv_team)).setOnItemClickListener(this);

        ((CompoundButton)findViewById(R.id.cb_cloth)).setOnCheckedChangeListener(this);
        ((CompoundButton)findViewById(R.id.cb_hat)).setOnCheckedChangeListener(this);
        ((CompoundButton)findViewById(R.id.cb_ball)).setOnCheckedChangeListener(this);

        String [] str_league = {"國家聯盟","美國聯盟"};
        ArrayAdapter<String> ad_league = new ArrayAdapter<>(this, android.R.layout.simple_spinner_dropdown_item,str_league);
        ((Spinner)findViewById(R.id.spn_league)).setAdapter(ad_league);

        String [] str_area = {"東區","中區","西區"};
        ArrayAdapter<String> ad_area = new ArrayAdapter<>(this, android.R.layout.simple_spinner_dropdown_item,str_area);
        ((Spinner)findViewById(R.id.spn_area)).setAdapter(ad_area);

        String [] str_seat = {"內野(USD$84)","外野(USD$37)","本壘後方(USD$232)"};
        ArrayAdapter<String> ad_seat = new ArrayAdapter<>(this, android.R.layout.simple_spinner_dropdown_item,str_seat);
        ((Spinner)findViewById(R.id.spn_seat)).setAdapter(ad_seat);

    }

    String str,str_league="國家聯盟",str_area="東區",str_team_="亞特蘭大勇士隊",str_seat="內野";
    int check_league=0,check_area=0,sum_seat=84,sum_buy=0,check=1,check_team=1;
    ArrayList<String> str_team = new ArrayList<>();

    @Override
    public void onClick(View view) {
        if(check!=check_team){
            Toast.makeText(this,"請選擇球隊",Toast.LENGTH_LONG).show();
        }else{
            new DatePickerDialog(this,this,2021,0,1).show();
            str="您選擇的時間為\n";
        }
    }

    @Override
    public void onDateSet(DatePicker datePicker, int i, int i1, int i2) {
        str+=String.format("%d年%d月%d日",i,i1+1,i2);
        new TimePickerDialog(this,this,12,00,true).show();
    }

    @Override
    public void onTimeSet(TimePicker timePicker, int i, int i1) {
        str+=String.format("%02d時%02d分\n",i,i1);
        AlertDialog.Builder bud= new AlertDialog.Builder(this);
                bud.setTitle("購票資訊");
                bud.setMessage(str + str_team_ + ",總金額USD$" + (sum_seat+sum_buy));
                bud.setNegativeButton("取消",this);
                bud.setPositiveButton("確定",this);
                bud.setCancelable(false);
                bud.show();
    }

    @Override
    public void onClick(DialogInterface dialogInterface, int i) {
        if(i==-1){
            Snackbar.make(findViewById(R.id.root),"購票成功",Snackbar.LENGTH_LONG).show();
        }else if(i==-2){
            Toast.makeText(this,"取消購票",Toast.LENGTH_LONG).show();
        }
    }

    @Override
    public void onItemSelected(AdapterView<?> adapterView, View view, int i, long l) {
        if(adapterView.getId()==R.id.spn_league){
            if(i==0){
                str_league="國家聯盟";
            }else if(i==1){
                str_league="美國聯盟";
            }
            check_league=i;
            check_area=0;
        }

        if(adapterView.getId()==R.id.spn_area){
            if(i==0){
                if(check_league==0){
                    check=1;
                }else{
                    check=4;
                }
                str_area="東區";
            }else if(i==1){
                if(check_league==0){
                    check=2;
                }else{
                    check=5;
                }
                str_area="中區";
            }else if(i==2){
                if(check_league==0){
                    check=3;
                }else{
                    check=6;
                }
                str_area="西區";
            }
            check_area=i;
        }

        if(check_league==0){
            if(check_area==0){
                str_team.clear();
                str_team.add("亞特蘭大勇士隊");
                str_team.add("邁阿密馬林魚隊");
                str_team.add("紐約大都會隊");
                str_team.add("費城費城人隊");
                str_team.add("華盛頓國民隊");
                check_team=1;
            }else if(check_area==1){
                str_team.clear();
                str_team.add("芝加哥小熊隊");
                str_team.add("辛辛那提紅人隊");
                str_team.add("密爾瓦基釀酒人隊");
                str_team.add("匹茲堡海盜隊");
                str_team.add("聖路易紅雀隊");
                check_team=2;
            }else if(check_area==2){
                str_team.clear();
                str_team.add("亞利桑那響尾蛇隊");
                str_team.add("科羅拉多落磯隊");
                str_team.add("洛杉磯道奇隊");
                str_team.add("聖地牙哥教士隊");
                str_team.add("舊金山巨人隊");
                check_team=3;
            }
        }else if(check_league==1){
            if(check_area==0){
                str_team.clear();
                str_team.add("巴爾的摩金鶯隊");
                str_team.add("波士頓紅襪隊");
                str_team.add("紐約洋基隊");
                str_team.add("坦帕灣光芒隊");
                str_team.add("多倫多藍鳥隊");
                check_team=4;
            }else if(check_area==1){
                str_team.clear();
                str_team.add("芝加哥白襪隊");
                str_team.add("克里夫蘭印地安人隊");
                str_team.add("底特律老虎隊");
                str_team.add("堪薩斯市皇家隊");
                str_team.add("明尼蘇達雙城隊");
                check_team=5;
            }else if(check_area==2){
                str_team.clear();
                str_team.add("休士頓太空人隊");
                str_team.add("洛杉磯天使隊");
                str_team.add("奧克蘭運動家隊");
                str_team.add("西雅圖水手隊");
                str_team.add("德州遊騎兵隊");
                check_team=6;
            }
        }
        ArrayAdapter<String> ad_team = new ArrayAdapter<>(this, android.R.layout.simple_list_item_1,str_team);
        ((ListView)findViewById(R.id.lsv_team)).setAdapter(ad_team);

        sum_seat=0;
        if(adapterView.getId()==R.id.spn_seat){
            if(i==0){
                str_seat = "內野";
                sum_seat += 84;
            }else if(i==1){
                str_seat = "外野";
                sum_seat += 37;
            }else if(i==2){
                str_seat = "本壘後方";
                sum_seat += 232;
           }
        }
        ((TextView)findViewById(R.id.txv_show)).setText(str_league+" "+str_area+" "+str_team_+" "+str_seat);
    }

    @Override
    public void onNothingSelected(AdapterView<?> adapterView) {

    }

    @Override
    public void onItemClick(AdapterView<?> adapterView, View view, int i, long l) {
        str_team_ = ((ListView)findViewById(R.id.lsv_team)).getItemAtPosition(i).toString();
        ((TextView)findViewById(R.id.txv_show)).setText(str_league+" "+str_area+" "+str_team_+" "+str_seat);
    }

    @Override
    public void onCheckedChanged(CompoundButton compoundButton, boolean b) {
        sum_buy = 0;
        if(((CompoundButton)findViewById(R.id.cb_cloth)).isChecked()){
            sum_buy += 50;
        }
        if(((CompoundButton)findViewById(R.id.cb_hat)).isChecked()){
            sum_buy += 29;
        }
        if(((CompoundButton)findViewById(R.id.cb_ball)).isChecked()){
            sum_buy += 42;
        }
    }
}