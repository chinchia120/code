package ncku.gm.p1006;

import androidx.appcompat.app.AppCompatActivity;

import android.graphics.Color;
import android.os.Bundle;
import android.view.View;
import android.widget.EditText;
import android.widget.TextView;

public class MainActivity extends AppCompatActivity {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);
    }

    public void show (View v){
        EditText edt_fn=findViewById(R.id.edt_first_name);
        EditText edt_n=findViewById(R.id.edt_name);
        EditText edt_pn=findViewById(R.id.edt_phone_number);

        String s1=edt_fn.getText().toString();
        String s2=edt_n.getText().toString();
        String s3=edt_pn.getText().toString();

        TextView txv_show=findViewById(R.id.txv_msg);
        txv_show.setText(s1+s2+"的電話是\n"+s3);
        //txv_show.setTextColor(Color.RED);
    }
}