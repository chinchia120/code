package ncku.gm.practice_6;

import androidx.appcompat.app.AppCompatActivity;

import android.app.AlertDialog;
import android.app.DatePickerDialog;
import android.app.TimePickerDialog;
import android.content.DialogInterface;
import android.os.Bundle;
import android.view.View;
import android.widget.Button;
import android.widget.DatePicker;
import android.widget.TextView;
import android.widget.TimePicker;
import android.widget.Toast;

import com.google.android.material.snackbar.Snackbar;

public class MainActivity extends AppCompatActivity
        implements View.OnClickListener , DatePickerDialog.OnDateSetListener ,
        TimePickerDialog.OnTimeSetListener , DialogInterface.OnClickListener {

    String str;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        ((Button)findViewById(R.id.btn_book)).setOnClickListener(this);


    }

    @Override
    public void onClick(View view) {
        new DatePickerDialog(this,this,2021,0,1).show();
        str="您選擇的時間為\n";
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
                bud.setMessage(str);
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
}