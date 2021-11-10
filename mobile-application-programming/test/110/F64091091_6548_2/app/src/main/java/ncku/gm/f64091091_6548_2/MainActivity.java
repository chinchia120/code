package ncku.gm.f64091091_6548_2;

import androidx.appcompat.app.AlertDialog;
import androidx.appcompat.app.AppCompatActivity;

import android.content.DialogInterface;
import android.os.Bundle;
import android.view.View;
import android.widget.AdapterView;
import android.widget.ArrayAdapter;
import android.widget.Button;
import android.widget.ListView;
import android.widget.Spinner;
import android.widget.TextView;

import com.google.android.material.snackbar.Snackbar;

import java.util.ArrayList;

public class MainActivity extends AppCompatActivity
        implements View.OnClickListener , View.OnLongClickListener ,
        AdapterView.OnItemClickListener , DialogInterface.OnClickListener , AdapterView.OnItemSelectedListener {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        ((Button)findViewById(R.id.btn_ear_minus)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_ear_plus)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_hand)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_hand)).setOnLongClickListener(this);

        ((ListView)findViewById(R.id.lsv_buy)).setOnItemClickListener(this);

        ((Spinner)findViewById(R.id.spn_disk)).setOnItemSelectedListener(this);

        String [] str_disk = {"0","1","2","3","4"};
        ArrayAdapter<String> ad_disk = new ArrayAdapter<>(this, android.R.layout.simple_spinner_dropdown_item,str_disk);
        ((Spinner)findViewById(R.id.spn_disk)).setAdapter(ad_disk);

        String [] str_buy = {"貨到付款","線上刷卡"};
        ArrayAdapter<String> ad_buy = new ArrayAdapter<>(this, android.R.layout.simple_list_item_1,str_buy);
        ((ListView)findViewById(R.id.lsv_buy)).setAdapter(ad_buy);
    }

    int ear=0,hand=0,disk=0,flag=0,sum=0;
    String str_buy;

    @Override
    public void onClick(View view) {
        if(view.getId()==R.id.btn_ear_minus){
            ear--;
            check();
            if(flag==1){
                ear++;
            }
        }else if(view.getId()==R.id.btn_ear_plus){
            ear++;
        }else if(view.getId()==R.id.btn_hand){
            hand++;
        }
        show();
    }

    @Override
    public boolean onLongClick(View view) {
        hand--;
        check();
        if(flag==1){
            hand++;
        }
        show();
        return true;
    }

    public void check(){
        flag=0;
        if(ear<0 || hand <0){
            flag=1;
        }
    }

    @Override
    public void onItemClick(AdapterView<?> adapterView, View view, int i, long l) {
        if(i==0){
            str_buy = "貨到付款";
        }else if(i==1){
            str_buy = "線上刷卡";
        }
        AlertDialog.Builder bud = new AlertDialog.Builder(this);
        bud.setTitle("購物訊息");
        bud.setMessage("付款方式為" + str_buy + ",總金額為$"+sum);
        bud.setCancelable(false);
        bud.setPositiveButton("確認",this);
        bud.setNegativeButton("取消",this);
        bud.show();
    }

    @Override
    public void onClick(DialogInterface dialogInterface, int i) {
        if(i==-1){
            Snackbar.make(findViewById(R.id.root),"總金額為$"+sum,Snackbar.LENGTH_LONG).show();
        }else if(i==-2){
            Snackbar.make(findViewById(R.id.root),"取消訂單",Snackbar.LENGTH_LONG).show();
        }
    }

    public void show(){
        ((TextView)findViewById(R.id.txv_ear)).setText(ear+"");
        ((TextView)findViewById(R.id.txv_hand)).setText(hand+"");
        ((TextView)findViewById(R.id.txv_disk)).setText(disk+"");
    }

    @Override
    public void onItemSelected(AdapterView<?> adapterView, View view, int i, long l) {
        disk=i;
        show();
    }

    @Override
    public void onNothingSelected(AdapterView<?> adapterView) {

    }
}