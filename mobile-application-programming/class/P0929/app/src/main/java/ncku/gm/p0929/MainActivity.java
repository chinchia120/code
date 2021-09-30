package ncku.gm.p0929;

import androidx.appcompat.app.AppCompatActivity;

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
    //int size=20;
    public void enlarge(View v){
        EditText edt_entername=(EditText) findViewById(R.id.edt_entername);
        String str=edt_entername.getText().toString();
        TextView txv_msg=(TextView) findViewById(R.id.txv_msg);
        txv_msg.setText(str +",您好!!");
        //size+=4;
        //txv_msg.setTextSize(size);
    }
}