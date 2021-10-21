package ncku.gm.p1020;

import androidx.appcompat.app.AppCompatActivity;

import android.os.Bundle;
import android.view.View;
import android.widget.CheckBox;
import android.widget.CompoundButton;
import android.widget.ImageView;
import android.widget.RadioGroup;
import android.widget.TextView;

public class MainActivity extends AppCompatActivity
        implements RadioGroup.OnCheckedChangeListener , CompoundButton.OnCheckedChangeListener {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        RadioGroup rg = findViewById(R.id.rg_ticket);
        rg.setOnCheckedChangeListener(this);
        CheckBox cb_h = findViewById(R.id.cb_ham);
        cb_h.setOnCheckedChangeListener(this);
        CheckBox cb_f = findViewById(R.id.cb_fra);
        cb_f.setOnCheckedChangeListener(this);
        CheckBox cb_c = findViewById(R.id.cb_cola);
        cb_c.setOnCheckedChangeListener(this);
        CheckBox cb_s = findViewById(R.id.cb_sp);
        cb_s.setOnCheckedChangeListener(this);
    }

    void show(){
        CheckBox cb_h = findViewById(R.id.cb_ham);
        CheckBox cb_f = findViewById(R.id.cb_fra);
        CheckBox cb_c = findViewById(R.id.cb_cola);
        CheckBox cb_s = findViewById(R.id.cb_sp);
        ImageView iv1 = findViewById(R.id.op1);
        ImageView iv2 = findViewById(R.id.op2);
        ImageView iv3 = findViewById(R.id.op3);
        ImageView iv4 = findViewById(R.id.op4);

        String str = "您的餐點是 : \n";

        if(cb_h.isChecked()){
            str+="漢堡\n";
            iv1.setVisibility(View.VISIBLE);
        }else{
            iv1.setVisibility(View.GONE);
        }

        if(cb_f.isChecked()) {
            str += "薯條\n";
            iv2.setVisibility(View.VISIBLE);
        }else{
            iv2.setVisibility(View.GONE);
        }

        if(cb_c.isChecked()){
            str+="可樂\n";
            iv3.setVisibility(View.VISIBLE);
        }else{
            iv3.setVisibility(View.GONE);
        }

        if(cb_s.isChecked()){
            str+="玉米濃湯\n";
            iv4.setVisibility(View.VISIBLE);
        }else{
            iv4.setVisibility(View.GONE);
        }

        TextView txv = findViewById(R.id.txv_msg);
        txv.setText(str);

        /*TextView txv = findViewById(R.id.txv_msg);
        RadioGroup rg = findViewById(R.id.rg_ticket);

        if(rg.getCheckedRadioButtonId() == R.id.rbtn_adult){
            txv.setText("買全票");
        }else if(rg.getCheckedRadioButtonId() == R.id.rbtn_child){
            txv.setText("買半票");
        }else{
            txv.setText("買敬老票");
        }*/
    }

    public void ok(View v){
        show();
    }

    @Override
    public void onCheckedChanged(RadioGroup radioGroup, int i) {
        show();
    }

    @Override
    public void onCheckedChanged(CompoundButton compoundButton, boolean b) {
        show();
    }
}