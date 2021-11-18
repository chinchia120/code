package ncku.gm.p1117;

import androidx.appcompat.app.AppCompatActivity;

import android.content.Intent;
import android.os.Bundle;
import android.view.View;
import android.widget.Button;
import android.widget.EditText;
import android.widget.TextView;

public class MainActivity2 extends AppCompatActivity implements View.OnClickListener {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main2);

        ((Button)findViewById(R.id.btn_call_act1)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_save)).setOnClickListener(this);

        Intent it2 = getIntent();
        int id = it2.getIntExtra("編號",0);
        String text = it2.getStringExtra("內容");

        ((TextView)findViewById(R.id.txv_id)).setText(id+"");
        ((TextView)findViewById(R.id.txv_text)).setText(text);
        ((EditText)findViewById(R.id.edt_enter)).setText(text);
    }

    @Override
    public void onClick(View view) {
        if(view.getId()==R.id.btn_call_act1){
            //startActivity(new Intent(this,MainActivity.class));
            finish();
        }else if(view.getId()==R.id.btn_save){
            Intent it3 = new Intent();
            String str_id = ((TextView)findViewById(R.id.txv_id)).getText().toString();
            int id = Integer.parseInt(str_id);
            it3.putExtra("編號",id);

            //String str_enter = ((EditText)findViewById(R.id.edt_enter)).getText().toString();
            it3.putExtra("內容",((EditText)findViewById(R.id.edt_enter)).getText().toString());

            setResult(RESULT_OK,it3);
            finish();
        }
    }
}