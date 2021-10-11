package ncku.gm.practice_2;

import androidx.appcompat.app.AppCompatActivity;

import android.os.Bundle;
import android.text.TextUtils;
import android.view.View;
import android.widget.EditText;
import android.widget.TextView;

import java.util.ArrayList;

public class MainActivity extends AppCompatActivity {
    TextView txv_2,txv_10,txv_16;
    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        txv_2 =(TextView) findViewById(R.id.txv_2);
        txv_10 =(TextView) findViewById(R.id.txv_10);
        txv_16 =(TextView) findViewById(R.id.txv_16);
    }

    public void click_0(View v){
        txv_2.setText(txv_2.getText().toString()+"0");
    }

    public void click_1(View v) {
        txv_2.setText(txv_2.getText().toString()+"1");
    }

    public void delete(View v){
        String str=txv_2.getText().toString();
        if(str.length()==0){
            txv_2.setText("");
        }else{
            str=str.substring(0,str.length()-1);
            txv_2.setText(str);
        }
    }

    public void reset(View v){
        txv_2.setText("");
        txv_10.setText("十進位 = ");
        txv_16.setText("十六進位 = ");
    }

    public void calculate(View v){
        String str=txv_2.getText().toString();
        if(TextUtils.isEmpty(str)){
            txv_2.setText("");
            txv_10.setText("十進位 = ");
            txv_16.setText("十六進位 = ");
        }else{
            int sum_10=Integer.valueOf(str,2);
            String str_16=Integer.toHexString(sum_10);
            txv_10.setText("十進位 = " + sum_10);
            txv_16.setText("十六進位 = " + str_16);
        }
    }
}