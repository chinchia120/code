package ncku.gm.practice_3;

import androidx.appcompat.app.AppCompatActivity;

import android.content.Context;
import android.os.Bundle;
import android.os.Vibrator;
import android.text.TextUtils;
import android.view.View;
import android.widget.Button;
import android.widget.EditText;
import android.widget.TextView;

public class MainActivity extends AppCompatActivity {
    TextView txv;
    EditText edt;
    Button btn_pl,btn_mi;
    Vibrator vb;
    int cnt=0;
    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        txv = findViewById(R.id.txv_num);
        edt =findViewById(R.id.edt_enter);
        btn_pl = findViewById(R.id.btn_plus);
        btn_mi = findViewById(R.id.btn_minus);

        btn_pl.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                txv.setText(++cnt + "");
                vibration(view);
            }
        });

        btn_pl.setOnLongClickListener(new View.OnLongClickListener() {
            @Override
            public boolean onLongClick(View view) {
                cnt+=2;
                txv.setText(cnt + "");
                vibration(view);
                return true;
            }
        });

        btn_mi.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                txv.setText(--cnt + "");
                vibration(view);
            }
        });

        btn_mi.setOnLongClickListener(new View.OnLongClickListener() {
            @Override
            public boolean onLongClick(View view) {
                cnt-=2;
                txv.setText(cnt + "");
                vibration(view);
                return true;
            }
        });


    }
    public void vibration(View v){
        vb=(Vibrator) getSystemService(Context.VIBRATOR_SERVICE);
        String str=edt.getText().toString();
        if(TextUtils.isEmpty(str)){
            return;
        }
        int n=Integer.parseInt(str);
        if(cnt!=0){
            if(Math.abs(cnt)%n == 0){
                vb.vibrate(3000);
            }
        }
    }


}