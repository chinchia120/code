package ncku.gm.f64091091_0108_1;

import androidx.appcompat.app.AppCompatActivity;

import android.app.DatePickerDialog;
import android.os.Bundle;
import android.view.View;
import android.widget.AdapterView;
import android.widget.ArrayAdapter;
import android.widget.Button;
import android.widget.CheckBox;
import android.widget.DatePicker;
import android.widget.EditText;
import android.widget.ImageView;
import android.widget.RadioButton;
import android.widget.Spinner;
import android.widget.TextView;

public class MainActivity extends AppCompatActivity
        implements View.OnClickListener , DatePickerDialog.OnDateSetListener ,
        AdapterView.OnItemSelectedListener {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        ((Button)findViewById(R.id.btn_date)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_ok)).setOnClickListener(this);

        ((Spinner)findViewById(R.id.spn_study)).setOnItemSelectedListener(this);

        String [] str_study = {"國中","高中","大學","研究所"};
        ArrayAdapter<String> ad_study = new ArrayAdapter<>(this, android.R.layout.simple_spinner_dropdown_item,str_study);
        ((Spinner)findViewById(R.id.spn_study)).setAdapter(ad_study);
    }

    String str_show , str_name , str_blood ,str_study;
    int old;

    @Override
    public void onClick(View view) {
        if(view.getId()==R.id.btn_date){
            new DatePickerDialog(this,this,2001,0,1).show();
        }else if(view.getId()==R.id.btn_ok){
            data_set();
        }
    }

    @Override
    public void onDateSet(DatePicker datePicker, int i, int i1, int i2) {
        ((TextView)findViewById(R.id.txv_show_date)).setText(String.format("%d/%d/%d",i,i1+1,i2));
        old = 2021 - i ;
    }

    public void data_set(){
        str_name = ((EditText)findViewById(R.id.edt_name)).getText().toString();
        char str_first_name = str_name.charAt(0);
        char str_last_name_1 = str_name.charAt(1);
        char str_last_name_2 = str_name.charAt(2);

        if(((RadioButton)findViewById(R.id.rbtn_a)).isChecked()){
            str_blood="A";
        }else if(((RadioButton)findViewById(R.id.rbtn_b)).isChecked()){
            str_blood="B";
        }else if(((RadioButton)findViewById(R.id.rbtn_o)).isChecked()){
            str_blood="O";
        }else if(((RadioButton)findViewById(R.id.rbtn_ab)).isChecked()){
            str_blood="AB";
        }

        if(((CheckBox)findViewById(R.id.cb_basketball)).isChecked()){
            ((ImageView)findViewById(R.id.imv_basketball)).setVisibility(View.VISIBLE);
        }else{
            ((ImageView)findViewById(R.id.imv_basketball)).setVisibility(View.GONE);
        }

        if(((CheckBox)findViewById(R.id.cb_basketball)).isChecked()){
            ((ImageView)findViewById(R.id.imv_baseball)).setVisibility(View.VISIBLE);
        }else{
            ((ImageView)findViewById(R.id.imv_baseball)).setVisibility(View.GONE);
        }

        if(((CheckBox)findViewById(R.id.cb_basketball)).isChecked()){
            ((ImageView)findViewById(R.id.imv_badminton)).setVisibility(View.VISIBLE);
        }else{
            ((ImageView)findViewById(R.id.imv_badminton)).setVisibility(View.GONE);
        }

        str_show = "履歷 :\n我姓" + str_first_name + "名" + str_last_name_1 + str_last_name_2 + ",今年" + old + "歲"
        + ",血型" + str_blood + "型," + str_study + "畢業,興趣如下\n" ;
        ((TextView)findViewById(R.id.txv_show)).setText(str_show);
    }

    @Override
    public void onItemSelected(AdapterView<?> adapterView, View view, int i, long l) {
        if(i==0){
            str_study = "國中";
        }else if(i==1){
            str_study="高中";
        }else if(i==2){
            str_study="大學";
        }else if(i==3){
            str_study="研究所";
        }
    }

    @Override
    public void onNothingSelected(AdapterView<?> adapterView) {

    }
}