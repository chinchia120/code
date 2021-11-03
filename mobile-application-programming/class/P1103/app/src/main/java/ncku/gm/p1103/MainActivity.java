package ncku.gm.p1103;

import androidx.appcompat.app.AlertDialog;
import androidx.appcompat.app.AppCompatActivity;

import android.app.DatePickerDialog;
import android.app.TimePickerDialog;
import android.content.DialogInterface;
import android.os.Bundle;
import android.view.Gravity;
import android.view.View;
import android.widget.AdapterView;
import android.widget.DatePicker;
import android.widget.ListView;
import android.widget.TimePicker;
import android.widget.Toast;

import com.google.android.material.snackbar.Snackbar;

public class MainActivity extends AppCompatActivity
        implements AdapterView.OnItemClickListener , DialogInterface.OnClickListener ,
        DatePickerDialog.OnDateSetListener , TimePickerDialog.OnTimeSetListener{

    Toast tos;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        ((ListView)findViewById(R.id.lsv_question)).setOnItemClickListener(this);

        tos = Toast.makeText(this,"",Toast.LENGTH_SHORT);
        tos.setGravity(Gravity.TOP | Gravity.LEFT,100,50);

        AlertDialog.Builder bdr = new AlertDialog.Builder(this);
        bdr.setIcon(android.R.drawable.presence_away);
        bdr.setTitle("歡迎");
        bdr.setMessage("示範");
        bdr.setCancelable(false);
        bdr.setNegativeButton("負向鈕",this);
        bdr.setPositiveButton("正向鈕",this);
        bdr.setNeutralButton("中性鈕",this);
        bdr.show();
    }

    @Override
    public void onItemClick(AdapterView<?> adapterView, View view, int i, long l) {
        if(i==0){
            tos.setText("球門");
            tos.show();
            }else if(i==1){
            tos.setText("虧");
            tos.show();
        }else if(i==2){
            tos.setText("傻瓜");
            tos.show();
        }else if(i==3){
            Snackbar.make(findViewById(R.id.root),"瀑布",Snackbar.LENGTH_SHORT).show();
        }else if(i==4){
            Snackbar.make(findViewById(R.id.root),"環保署",Snackbar.LENGTH_SHORT).show();
        }else if(i==5){
            Snackbar.make(findViewById(R.id.root),"偷笑",Snackbar.LENGTH_SHORT).show();
        }else  if(i==6){
            DatePickerDialog dpd = new DatePickerDialog(this,this,2021,10,3);
            dpd.show();
        }else if(i==7){
            new TimePickerDialog(this,this,15,55,true).show();
        }
    }

    @Override
    public void onClick(DialogInterface dialogInterface, int i) {
        if(i==DialogInterface.BUTTON_POSITIVE){
            Snackbar.make(findViewById(R.id.root),"正",Snackbar.LENGTH_SHORT).show();
        }else if(i==DialogInterface.BUTTON_NEGATIVE){
            Snackbar.make(findViewById(R.id.root),"負",Snackbar.LENGTH_SHORT).show();
        }else if(i==DialogInterface.BUTTON_NEUTRAL){
            Snackbar.make(findViewById(R.id.root),"中",Snackbar.LENGTH_SHORT).show();
        }
    }

    @Override
    public void onDateSet(DatePicker datePicker, int i, int i1, int i2) {
        String str = "日期 : "+i+"/"+(i1+1)+"/"+i2;
        Snackbar.make(findViewById(R.id.root),str,Snackbar.LENGTH_LONG).show();
    }

    @Override
    public void onTimeSet(TimePicker timePicker, int i, int i1) {
        String str = "時間 : "+String.format("%02d:%02d",i,i1);
        Snackbar.make(findViewById(R.id.root),str,Snackbar.LENGTH_LONG).show();
    }
}