package ncku.gm.f64091091_0108_2;

import androidx.appcompat.app.AppCompatActivity;

import android.app.AlertDialog;
import android.app.DatePickerDialog;
import android.content.DialogInterface;
import android.os.Bundle;
import android.view.View;
import android.widget.AdapterView;
import android.widget.ArrayAdapter;
import android.widget.Button;
import android.widget.ListView;
import android.widget.Spinner;
import android.widget.TextView;
import android.widget.Toast;

import java.util.ArrayList;

public class MainActivity extends AppCompatActivity
        implements View.OnClickListener , View.OnLongClickListener  ,
        AdapterView.OnItemClickListener , DialogInterface.OnClickListener {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        ((Button)findViewById(R.id.btn_num_in)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_num_in)).setOnLongClickListener(this);
        ((Button)findViewById(R.id.btn_minus1)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_plus1)).setOnClickListener(this);

        ((ListView)findViewById(R.id.lsv_buy)).setOnItemClickListener(this);

        tos=Toast.makeText(this,"不可大於6或小於0",Toast.LENGTH_LONG);

        String[] str_hat = {"0","1","2","3","4","5","6"};
        ArrayAdapter<String> ad_hat = new ArrayAdapter<>(this, android.R.layout.simple_spinner_dropdown_item,str_hat);
        ((Spinner)findViewById(R.id.spn_hat)).setAdapter(ad_hat);
    }

    int in=0,out=0,sum=0,flag=0;
    String str_show,str_buy;
    Toast tos;

    @Override
    public void onClick(View view) {
        if(view.getId()==R.id.btn_num_in){
            in++;
            check();
            if(flag==1){
                in--;
            }
        }else if(view.getId()==R.id.btn_minus1){
            out--;
            check();
            if(flag==1){
                out++;
            }
        }else if(view.getId()==R.id.btn_plus1){
            out++;
            check();
            if(flag==1){
                out--;
            }
        }
        show();
        ArrayList<String> str_hat = new ArrayList<>();
        str_hat.add("0");
        str_hat.add("1");
        str_hat.add("2");
        str_hat.add("3");
        str_hat.add("4");
        str_hat.add("5");
        str_hat.add("6");
        for(int i=6;i>(in+out);i--){
            str_hat.remove(i);
        }
        ArrayAdapter<String> ad_hat = new ArrayAdapter<>(this, android.R.layout.simple_spinner_dropdown_item,str_hat);
        ((Spinner)findViewById(R.id.spn_hat)).setAdapter(ad_hat);
    }

    @Override
    public boolean onLongClick(View view) {
        in--;
        check();
        if(flag==1){
            in++;
        }
        show();
        return true;
    }

    @Override
    public void onItemClick(AdapterView<?> adapterView, View view, int i, long l) {
        if(i==0){
            str_buy="現場付款";
        }else{
            str_buy="線上刷卡";
        }
        data_set();
        AlertDialog.Builder bud = new AlertDialog.Builder(this);
        bud.setTitle("購票資訊");
        bud.setMessage(str_show);
        bud.setNegativeButton("取消",this);
        bud.setPositiveButton("確認",this);
        bud.setCancelable(false);
        bud.show();
    }

    @Override
    public void onClick(DialogInterface dialogInterface, int i) {
        if(i==-1){
            ((TextView)findViewById(R.id.txv_show)).setText("總金額為$"+sum);
        }else if(i==-2){
            ((TextView)findViewById(R.id.txv_show)).setText("取消訂單");
        }
    }

    public void data_set(){
        str_show = "付款方式為" + str_buy + ",總金額為$" + sum ;
    }

    public void check(){
        flag=0;
        if((in+out>6) || in>6 || out>6 || in+out<0 || in<0 ||out<0){
            tos.show();
            flag=1;
        }
    }

    public void show(){
        ((TextView)findViewById(R.id.txv_show_in)).setText(in+"");
        ((TextView)findViewById(R.id.txv_show_out)).setText(out+"");
    }
}