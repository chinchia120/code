package ncku.gm.f64091091_6548_1;

import androidx.appcompat.app.AppCompatActivity;

import android.app.DatePickerDialog;
import android.os.Bundle;
import android.view.View;
import android.widget.AdapterView;
import android.widget.ArrayAdapter;
import android.widget.Button;
import android.widget.CompoundButton;
import android.widget.DatePicker;
import android.widget.EditText;
import android.widget.RadioGroup;
import android.widget.Spinner;
import android.widget.TextView;

import com.google.android.material.datepicker.OnSelectionChangedListener;

import java.util.ArrayList;

public class MainActivity extends AppCompatActivity
        implements View.OnClickListener , DatePickerDialog.OnDateSetListener ,
        RadioGroup.OnCheckedChangeListener  , CompoundButton.OnCheckedChangeListener  ,
        AdapterView.OnItemSelectedListener {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        ((Button)findViewById(R.id.btn_date)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_ok)).setOnClickListener(this);

        ((RadioGroup)findViewById(R.id.rg_get)).setOnCheckedChangeListener(this);

        ((CompoundButton)findViewById(R.id.cb_bear)).setOnCheckedChangeListener(this);
        ((CompoundButton)findViewById(R.id.cb_trip)).setOnCheckedChangeListener(this);
        ((CompoundButton)findViewById(R.id.cb_farm)).setOnCheckedChangeListener(this);
        ((CompoundButton)findViewById(R.id.cb_sport)).setOnCheckedChangeListener(this);
        ((CompoundButton)findViewById(R.id.cb_food)).setOnCheckedChangeListener(this);

        ((Spinner)findViewById(R.id.spn_location)).setOnItemSelectedListener(this);

    }

    String str_area,str_id,str_date,str_get,str_location_,str_test,str_show;
    int cnt_date=0,cnt_location=0;

    @Override
    public void onClick(View view) {
        if(view.getId()==R.id.btn_date){
            new DatePickerDialog(this,this,2021,10,10).show();
        }else if(view.getId()==R.id.btn_ok){
            show();
        }
    }

    @Override
    public void onDateSet(DatePicker datePicker, int i, int i1, int i2) {
        if(i1==10){
            str_date = 10-i2 + "";
            cnt_date = 10-i2;
        }else if(i1==9){
            str_date = 41-i2 + "";
            cnt_date = 41-i2;
        }
        ((TextView)findViewById(R.id.txv_date)).setText(String.format("%d/%d/%d",i-1911,i1+1,i2));
    }

    @Override
    public void onCheckedChanged(RadioGroup radioGroup, int i) {
        ArrayList<String> str_location = new ArrayList<>();

        if(radioGroup.getCheckedRadioButtonId()==R.id.rbtn_paper){
            str_location.clear();
            str_location.add("超商");
            str_location.add("郵局");
            str_get = "紙本領取";
        }else if(radioGroup.getCheckedRadioButtonId()==R.id.rbtn_digital){
            str_location.clear();
            str_location.add("信用卡");
            str_location.add("行動支付");
            str_location.add("電子票證");
            str_get = "數位綁定";
        }
        ArrayAdapter<String> ad_location = new ArrayAdapter<>(this, android.R.layout.simple_spinner_dropdown_item,str_location);
        ((Spinner)findViewById(R.id.spn_location)).setAdapter(ad_location);
    }

    public void show(){
        str_test = ((EditText)findViewById(R.id.edt_id)).getText().toString();
        if(str_test.charAt(0)=='A'){
            str_area = "您為台北市民";
        }else{
            str_area = "您不為台北市民";
        }
        str_id = str_test.substring(7,10);
        str_location_ = ((Spinner)findViewById(R.id.spn_location)).getSelectedItem().toString();

        str_show = str_area + ",身分證末三碼為" +  str_id + ",提早" + str_date + "日,以" + str_location_ + str_get + "預定振興五倍劵" ;
        ((TextView)findViewById(R.id.txv_show)).setText(str_show);
    }

    @Override
    public void onCheckedChanged(CompoundButton compoundButton, boolean b) {
        int flag=0;
        if(((CompoundButton)findViewById(R.id.cb_bear)).isChecked()){
            if(str_test.charAt(0)=='A'){
                ((TextView)findViewById(R.id.txv_bear)).setText("中獎");
                flag=1;
            }else{
                ((TextView)findViewById(R.id.txv_bear)).setText("未中獎");
            }
        }else{
            ((TextView)findViewById(R.id.txv_bear)).setText("未勾選");
        }

        if(((CompoundButton)findViewById(R.id.cb_trip)).isChecked()){
            if(str_id == "168" || str_id == "076"){
                ((TextView)findViewById(R.id.txv_trip)).setText("中獎");
                flag=1;
            }else{
                ((TextView)findViewById(R.id.txv_trip)).setText("未中獎");
            }
        }else{
            ((TextView)findViewById(R.id.txv_trip)).setText("未勾選");
        }

        if(((CompoundButton)findViewById(R.id.cb_farm)).isChecked()){
            if(cnt_date>10){
                ((TextView)findViewById(R.id.txv_farm)).setText("中獎");
                flag=1;
            }else{
                ((TextView)findViewById(R.id.txv_farm)).setText("未中獎");
            }
        }else{
            ((TextView)findViewById(R.id.txv_farm)).setText("未勾選");
        }

        if(((CompoundButton)findViewById(R.id.cb_sport)).isChecked()){
            if(cnt_location==2){
                ((TextView)findViewById(R.id.txv_sport)).setText("中獎");
                flag=1;
            }else{
                ((TextView)findViewById(R.id.txv_sport)).setText("未中獎");
            }
        }else{
            ((TextView)findViewById(R.id.txv_sport)).setText("未勾選");
        }

        if(((CompoundButton)findViewById(R.id.cb_food)).isChecked()){
            if(flag==0){
                ((TextView)findViewById(R.id.txv_food)).setText("中獎");
            }else{
                ((TextView)findViewById(R.id.txv_food)).setText("未中獎");
            }
        }else{
            ((TextView)findViewById(R.id.txv_food)).setText("未勾選");
        }
    }

    @Override
    public void onItemSelected(AdapterView<?> adapterView, View view, int i, long l) {
        cnt_location = i;
        show();
    }

    @Override
    public void onNothingSelected(AdapterView<?> adapterView) {

    }
}