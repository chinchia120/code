package ncku.gm.f64091091_0109_2;

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
import android.widget.DatePicker;
import android.widget.ListView;
import android.widget.Spinner;
import android.widget.TextView;
import android.widget.TimePicker;
import android.widget.Toast;

import java.util.ArrayList;

public class MainActivity extends AppCompatActivity
        implements View.OnClickListener  , DatePickerDialog.OnDateSetListener ,
        TimePickerDialog.OnTimeSetListener , DialogInterface.OnClickListener ,
        AdapterView.OnItemClickListener , AdapterView.OnItemSelectedListener {

    Toast tos;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        ((Button)findViewById(R.id.btn_date)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_time)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_minus1)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_plus1)).setOnClickListener(this);

        ((ListView)findViewById(R.id.lsv_buy)).setOnItemClickListener(this);

        ((Spinner)findViewById(R.id.spn_ticket_sis)).setOnItemSelectedListener(this);

        tos=Toast.makeText(this,"不可大於6或小於0",Toast.LENGTH_LONG);

        ArrayList<String> ticket_sis = new ArrayList<>();
        ticket_sis.add("0");
        ticket_sis.add("1");
        ticket_sis.add("2");
        ticket_sis.add("3");
        ticket_sis.add("4");
        ticket_sis.add("5");
        ticket_sis.add("6");
        ArrayAdapter<String> ad_ticket_sis = new ArrayAdapter<>(this,android.R.layout.simple_spinner_dropdown_item,ticket_sis);
        ((Spinner)findViewById(R.id.spn_ticket_sis)).setAdapter(ad_ticket_sis);
    }

    int ticket_bro=0,sum=0;
    String str_show,str_date,str_time,str_bro,str_sis,str_buy;

    @Override
    public void onClick(View view) {
        if(view.getId()==R.id.btn_date){
            new DatePickerDialog(this,this,2020,10,11).show();
        }
        if(view.getId()==R.id.btn_time){
            new TimePickerDialog(this,this,12,00,true).show();
        }

        if(view.getId()==R.id.btn_minus1){
            ticket_bro--;
        }
        if(view.getId()==R.id.btn_plus1){
            ticket_bro++;
        }
        if(ticket_bro>6 || ticket_bro<0){
            tos.show();
            if(ticket_bro>6){
                ticket_bro=6;
            }else{
                ticket_bro=0;
            }
        }
        ((TextView)findViewById(R.id.txv_ticket_bro)).setText(ticket_bro+"");
        ArrayList<String> ticket_sis = new ArrayList<>();
        ticket_sis.add("0");
        ticket_sis.add("1");
        ticket_sis.add("2");
        ticket_sis.add("3");
        ticket_sis.add("4");
        ticket_sis.add("5");
        ticket_sis.add("6");
        for(int i=6;i>ticket_bro;i--){
            ticket_sis.remove(i);
        }
        ArrayAdapter<String> ad_ticket_sis = new ArrayAdapter<>(this,android.R.layout.simple_spinner_dropdown_item,ticket_sis);
        ((Spinner)findViewById(R.id.spn_ticket_sis)).setAdapter(ad_ticket_sis);
    }

    @Override
    public void onDateSet(DatePicker datePicker, int i, int i1, int i2) {
        str_date = (i1+1) +"/" +i2;
        ((TextView)findViewById(R.id.txv_show_date)).setText(i+"年"+(i1+1)+"月"+i2+"日");
    }

    @Override
    public void onTimeSet(TimePicker timePicker, int i, int i1) {
        if(i>6 && i<19){
            str_time = "日間";
        }else{
            str_time = "夜間";
        }
        ((TextView)findViewById(R.id.txv_show_time)).setText(String.format("%02d : %02d",i,i1));
    }

    @Override
    public void onClick(DialogInterface dialogInterface, int i) {
        if(i==-1){
            ((TextView)findViewById(R.id.txv_show)).setText("總金額為$"+sum);
        }else if(i==-2){
            ((TextView)findViewById(R.id.txv_show)).setText("取消訂單");
        }
    }

    @Override
    public void onItemClick(AdapterView<?> adapterView, View view, int i, long l) {
        data_set();
        if(i==0){
            str_buy = "線上刷卡,\n總金額為$";
        }else{
            str_buy = "現場付款,\n總金額為$";
        }

        AlertDialog.Builder bud = new AlertDialog.Builder(this);
        bud.setMessage(str_show+str_buy+sum);
        bud.setNegativeButton("取消",this);
        bud.setPositiveButton("確認",this);
        bud.setCancelable(false);
        bud.show();
    }

    public void data_set(){
        str_bro = ((TextView)findViewById(R.id.txv_ticket_bro)).getText().toString();
        str_show = str_date + str_time + "簽名會\n" + "炭治郎*" + str_bro + "," + "彌豆子*" + str_sis + "付款方式為";
    }

    @Override
    public void onItemSelected(AdapterView<?> adapterView, View view, int i, long l) {
        str_sis = i + "";
        if(ticket_bro+i>6){
            tos.show();
            ((Spinner)findViewById(R.id.spn_ticket_sis)).setSelection(0);
        }
    }

    @Override
    public void onNothingSelected(AdapterView<?> adapterView) {

    }
}