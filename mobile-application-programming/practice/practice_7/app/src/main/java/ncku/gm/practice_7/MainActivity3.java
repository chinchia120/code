package ncku.gm.practice_7;

import androidx.appcompat.app.AlertDialog;
import androidx.appcompat.app.AppCompatActivity;

import android.app.TimePickerDialog;
import android.content.DialogInterface;
import android.content.Intent;
import android.os.Bundle;
import android.view.View;
import android.widget.ArrayAdapter;
import android.widget.Button;
import android.widget.EditText;
import android.widget.Spinner;
import android.widget.TextView;
import android.widget.TimePicker;

import com.google.android.material.snackbar.Snackbar;

import java.util.ArrayList;

public class MainActivity3 extends AppCompatActivity
        implements View.OnClickListener , DialogInterface.OnClickListener ,
        TimePickerDialog.OnTimeSetListener {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main3);

        ((Button)findViewById(R.id.btn_sure)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_time)).setOnClickListener(this);

        ArrayList<String> str_res = (ArrayList<String>) getIntent().getBundleExtra("餐廳").getSerializable("ArrayList");
        ArrayAdapter<String> ad_res = new ArrayAdapter<>(this, android.R.layout.simple_spinner_dropdown_item,str_res);
        ((Spinner)findViewById(R.id.spn_restaurant)).setAdapter(ad_res);
    }

    String str_phone,str_res="",str_time;

    @Override
    public void onClick(View view) {
        if(view.getId()==R.id.btn_sure){
            str_phone = ((EditText)findViewById(R.id.edt_phone)).getText().toString();
            str_res = ((Spinner)findViewById(R.id.spn_restaurant)).getSelectedItem().toString();
            AlertDialog.Builder bud = new AlertDialog.Builder(this);
            bud.setTitle("訂餐資訊");
            bud.setMessage("訂購 : "+ str_res +"\n連絡電話 : "+str_phone+"\n取餐時間 : "+str_time);
            bud.setNegativeButton("取消",this);
            bud.setPositiveButton("確定",this);
            bud.setCancelable(false);
            bud.show();
        }else if(view.getId()==R.id.btn_time){
            new TimePickerDialog(this,this,12,00,true).show();
        }
    }


    @Override
    public void onClick(DialogInterface dialogInterface, int i) {
        setResult(-1,new Intent().putExtra("確認",i));
        finish();
    }

    @Override
    public void onTimeSet(TimePicker timePicker, int i, int i1) {
        str_time = String.format("%02d:%02d",i,i1);
        ((TextView)findViewById(R.id.txv_time)).setText(str_time);
    }
}