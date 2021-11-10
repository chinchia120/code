package ncku.gm.f64091091_0109_1;

import androidx.appcompat.app.AppCompatActivity;

import android.content.Context;
import android.os.Bundle;
import android.os.Vibrator;
import android.view.View;
import android.widget.AdapterView;
import android.widget.ArrayAdapter;
import android.widget.Button;
import android.widget.CheckBox;
import android.widget.CompoundButton;
import android.widget.ImageView;
import android.widget.RadioButton;
import android.widget.RadioGroup;
import android.widget.Spinner;
import android.widget.TextView;

public class MainActivity extends AppCompatActivity
        implements View.OnClickListener , View.OnLongClickListener , AdapterView.OnItemSelectedListener ,
        RadioGroup.OnCheckedChangeListener , CompoundButton.OnCheckedChangeListener {

    Vibrator vb;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        ((Button)findViewById(R.id.btn_1)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_5)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_10)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_sick)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_1)).setOnLongClickListener(this);
        ((Button)findViewById(R.id.btn_5)).setOnLongClickListener(this);
        ((Button)findViewById(R.id.btn_10)).setOnLongClickListener(this);

        ((Spinner)findViewById(R.id.spn_name)).setOnItemSelectedListener(this);

        ((RadioGroup)findViewById(R.id.rg_blood)).setOnCheckedChangeListener(this);

        ((CompoundButton)findViewById(R.id.cb_china)).setOnCheckedChangeListener(this);
        ((CompoundButton)findViewById(R.id.cb_usa)).setOnCheckedChangeListener(this);
        ((CompoundButton)findViewById(R.id.cb_phi)).setOnCheckedChangeListener(this);
        ((CompoundButton)findViewById(R.id.cb_korea)).setOnCheckedChangeListener(this);

        String [] name = {"習大","川普","老杜","恩恩"};
        ArrayAdapter<String> ad_name = new ArrayAdapter<>(this, android.R.layout.simple_spinner_dropdown_item,name);
        ((Spinner)findViewById(R.id.spn_name)).setAdapter(ad_name);
    }

    int old = 0,cnt_china=0,cnt_usa=0,cnt_phi=0,cnt_korea=0;
    String str_show,str_name,str_old,str_blood="A",str_visit=" ";

    @Override
    public void onClick(View view) {
        if(view.getId()==R.id.btn_10){
            old+=10;
        }else if(view.getId()==R.id.btn_5){
            old+=5;
        }else if(view.getId()==R.id.btn_1){
            old+=1;
        }
        if(old<0){
            old=0;
        }else if(old>99){
            old=99;
        }
        ((TextView)findViewById(R.id.txv_showold)).setText(old+"");
        show();

        vb =(Vibrator) getSystemService(Context.VIBRATOR_SERVICE);
        if(view.getId()==R.id.btn_sick){
            if(((CompoundButton)findViewById(R.id.cb_china)).isChecked()){
                cnt_china++;
            }
            if(((CompoundButton)findViewById(R.id.cb_usa)).isChecked()){
                cnt_usa++;
            }
            if(((CompoundButton)findViewById(R.id.cb_phi)).isChecked()){
                cnt_phi++;
            }
            if(((CompoundButton)findViewById(R.id.cb_korea)).isChecked()){
                cnt_korea++;
            }
            if(cnt_china>0){
                ((TextView)findViewById(R.id.txv_cnt_china)).setText(cnt_china+"");
                ((ImageView)findViewById(R.id.imv_china)).setVisibility(View.VISIBLE);
            }
            if(cnt_usa>0){
                ((TextView)findViewById(R.id.txv_cnt_usa)).setText(cnt_usa+"");
                ((ImageView)findViewById(R.id.imv_usa)).setVisibility(View.VISIBLE);
            }
            if(cnt_phi>0){
                ((TextView)findViewById(R.id.txv_cnt_phi)).setText(cnt_phi+"");
                ((ImageView)findViewById(R.id.imv_phi)).setVisibility(View.VISIBLE);
            }
            if(cnt_korea>0){
                ((TextView)findViewById(R.id.txv_cnt_korea)).setText(cnt_korea+"");
                ((ImageView)findViewById(R.id.imv_korea)).setVisibility(View.VISIBLE);
            }
            if(cnt_china>2 || cnt_korea>2 || cnt_phi>2 || cnt_usa>2 ){
                vb.vibrate(3000);
            }
        }
    }

    @Override
    public boolean onLongClick(View view) {
        if(view.getId()==R.id.btn_10){
            old-=10;
        }else if(view.getId()==R.id.btn_5){
            old-=5;
        }else if(view.getId()==R.id.btn_1){
            old-=1;
        }
        if(old<0){
            old=0;
        }else if(old>99){
            old=99;
        }
        ((TextView)findViewById(R.id.txv_showold)).setText(old+"");
        show();
        return true;
    }

    @Override
    public void onItemSelected(AdapterView<?> adapterView, View view, int i, long l) {
        show();
    }

    @Override
    public void onNothingSelected(AdapterView<?> adapterView) {

    }

    @Override
    public void onCheckedChanged(RadioGroup radioGroup, int i) {
        if(radioGroup.getCheckedRadioButtonId()==R.id.rbtn_A){
            str_blood="A";
        }else if(radioGroup.getCheckedRadioButtonId()==R.id.rbtn_B){
            str_blood="B";
        }else if(radioGroup.getCheckedRadioButtonId()==R.id.rbtn_O){
            str_blood="O";
        }else if(radioGroup.getCheckedRadioButtonId()==R.id.rbtn_AB){
            str_blood="AB";
        }
        show();
    }

    @Override
    public void onCheckedChanged(CompoundButton compoundButton, boolean b) {
        str_visit=" \n";
        if(((CompoundButton)findViewById(R.id.cb_china)).isChecked()){
            str_visit+="中國 ";
        }
        if(((CompoundButton)findViewById(R.id.cb_usa)).isChecked()){
            str_visit+="美國 ";
        }
        if(((CompoundButton)findViewById(R.id.cb_phi)).isChecked()){
            str_visit+="菲律賓 ";
        }
        if(((CompoundButton)findViewById(R.id.cb_korea)).isChecked()){
            str_visit+="北韓 ";
        }
        show();
    }
    
    public void show(){
        str_name = ((Spinner)findViewById(R.id.spn_name)).getSelectedItem().toString();
        str_old = ((TextView)findViewById(R.id.txv_showold)).getText().toString();
        str_show="我是"+str_name+","+"今年"+str_old+"歲,"+"血型"+str_blood+"型,"+"曾造訪"+str_visit;
        ((TextView)findViewById(R.id.txv_msg)).setText(str_show+"");
    }
}