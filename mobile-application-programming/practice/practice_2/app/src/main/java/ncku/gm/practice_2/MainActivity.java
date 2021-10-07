package ncku.gm.practice_2;

import androidx.appcompat.app.AppCompatActivity;

import android.os.Bundle;
import android.view.View;
import android.widget.EditText;
import android.widget.TextView;

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

    }

    public void reset(View v){
        txv_2.setText("");
    }

    public void calculate(View v){
        int num=Integer.parseInt(txv_2.toString());

        txv_10.setText(num);
    }
}