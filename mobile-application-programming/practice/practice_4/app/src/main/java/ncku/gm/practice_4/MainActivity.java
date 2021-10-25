package ncku.gm.practice_4;

import androidx.appcompat.app.AppCompatActivity;

import android.os.Bundle;
import android.view.View;
import android.widget.Button;
import android.widget.CheckBox;
import android.widget.CompoundButton;
import android.widget.ImageView;
import android.widget.RadioButton;
import android.widget.RadioGroup;
import android.widget.TextView;

public class MainActivity extends AppCompatActivity {

    TextView txv;
    RadioButton rb_h,rb_f,rb_rt,rb_c,rb_hs,rb_ns,rb_li,rb_ni;
    CheckBox cb_f,cb_i;
    Button btn;
    RadioGroup rg_m,rg_d,rg_s,rg_i;
    ImageView imv1,imv2,imv3,imv4,imv5,imv6;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        txv = findViewById(R.id.txv_msg);

        rb_h = findViewById(R.id.rbtn_hamburger);
        rb_f = findViewById(R.id.rbtn_chicken);
        rb_rt = findViewById(R.id.rbtn_redtea);
        rb_c = findViewById(R.id.rbtn_coffee);
        rb_hs = findViewById(R.id.rbtn_halfsugar);
        rb_ns = findViewById(R.id.rbtn_nosugar);
        rb_li = findViewById(R.id.rbtn_lessice);
        rb_ni = findViewById(R.id.rbtn_noice);

        cb_f = findViewById(R.id.cb_frenchfries);
        cb_i = findViewById(R.id.cb_ice);

        btn = findViewById(R.id.btn_ok);

        rg_m = findViewById(R.id.rg_main);
        rg_d = findViewById(R.id.rg_drink);
        rg_s = findViewById(R.id.rg_sugar);
        rg_i = findViewById(R.id.rg_ice);

        imv1 = findViewById(R.id.imv_hamburger);
        imv2 = findViewById(R.id.imv_chicken);
        imv3 = findViewById(R.id.imv_redtea);
        imv4 = findViewById(R.id.imv_coffee);
        imv5 = findViewById(R.id.imv_frenchfeies);
        imv6 = findViewById(R.id.imv_ice);

        rg_m.setOnCheckedChangeListener(new RadioGroup.OnCheckedChangeListener() {
            @Override
            public void onCheckedChanged(RadioGroup radioGroup, int i) {
                show();
            }
        });

        rg_d.setOnCheckedChangeListener(new RadioGroup.OnCheckedChangeListener() {
            @Override
            public void onCheckedChanged(RadioGroup radioGroup, int i) {
                show();
            }
        });

        rg_s.setOnCheckedChangeListener(new RadioGroup.OnCheckedChangeListener() {
            @Override
            public void onCheckedChanged(RadioGroup radioGroup, int i) {
                show();
            }
        });

        rg_i.setOnCheckedChangeListener(new RadioGroup.OnCheckedChangeListener() {
            @Override
            public void onCheckedChanged(RadioGroup radioGroup, int i) {
                show();
            }
        });

        cb_f.setOnCheckedChangeListener(new CompoundButton.OnCheckedChangeListener() {
            @Override
            public void onCheckedChanged(CompoundButton compoundButton, boolean b) {
                show();
            }
        });

        cb_i.setOnCheckedChangeListener(new CompoundButton.OnCheckedChangeListener() {
            @Override
            public void onCheckedChanged(CompoundButton compoundButton, boolean b) {
                show();
            }
        });

        btn.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                show();
            }
        });
    }

    public void show(){

        String str = "點餐為 : \n";
        int sum=0;

        if(rb_h.isChecked()){
            str+="漢堡";
            sum+=40;
            imv1.setVisibility(View.VISIBLE);
            imv2.setVisibility(View.GONE);
        }else{
            str+="炸雞";
            sum+=50;
            imv1.setVisibility(View.GONE);
            imv2.setVisibility(View.VISIBLE);
        }

        if(rb_rt.isChecked()){
            str+="搭配紅茶";
            sum+=20;
            imv3.setVisibility(View.VISIBLE);
            imv4.setVisibility(View.GONE);

            if(rb_hs.isChecked()){
                str+="(半糖、";
            }else{
                str+="(無糖、";
            }

            if(rb_li.isChecked()){
                str+="少冰)，";
            }else{
                str+="去冰)，";
            }
        }else if(rb_c.isChecked()){
            str+="搭配咖啡";
            sum+=30;
            imv3.setVisibility(View.GONE);
            imv4.setVisibility(View.VISIBLE);

            if(rb_hs.isChecked()){
                str+="(半糖、";
            }else{
                str+="(無糖、";
            }

            if(rb_li.isChecked()){
                str+="少冰)，";
            }else{
                str+="去冰)，";
            }
        }

        if(cb_f.isChecked()){
            if(cb_i.isChecked()){
                str+="加點薯條和聖代，";
                sum+=50;
                imv5.setVisibility(View.VISIBLE);
                imv6.setVisibility(View.VISIBLE);
            }else {
                str+="加點薯條，";
                sum+=20;
                imv5.setVisibility(View.VISIBLE);
                imv6.setVisibility(View.GONE);
            }
        }else if(cb_i.isChecked()){
            str+="加點聖代，";
            sum+=30;
            imv5.setVisibility(View.GONE);
            imv6.setVisibility(View.VISIBLE);
        }

        str+="總共";
        str+=Integer.toString(sum);
        str+="元。";
        txv.setText(str);

    }
}