package ncku.gm.practice_1;

import androidx.appcompat.app.AppCompatActivity;

import android.graphics.Typeface;
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

    public void sayhello(View v){
        EditText name=(EditText) findViewById(R.id.edt_entername);
        TextView txv=(TextView) findViewById(R.id.txv_msg);
        txv.setText(name.getText().toString()+",您好!!");
    }

    int size=34;
    public void enlarge(View v){
        TextView txv=(TextView) findViewById(R.id.txv_msg);
        size+=3;
        txv.setTextSize(size);
    }

    public void shrink(View v){
        TextView txv=(TextView) findViewById(R.id.txv_msg);
        size-=3;
        txv.setTextSize(size);
    }

    public void red(View v){
        TextView txv=(TextView) findViewById(R.id.txv_msg);
        txv.setTextColor(this.getResources().getColor(R.color.red));
    }

    public void blue(View v){
        TextView txv=(TextView) findViewById(R.id.txv_msg);
        txv.setTextColor(this.getResources().getColor(R.color.blue));
    }

    public void bold(View v){
        TextView txv=(TextView) findViewById(R.id.txv_msg);
        Typeface font=Typeface.create(Typeface.SANS_SERIF,Typeface.BOLD);
        txv.setTypeface(font);
    }

    public void italic(View v){
        TextView txv=(TextView) findViewById(R.id.txv_msg);
        Typeface font=Typeface.create(Typeface.SANS_SERIF,Typeface.ITALIC);
        txv.setTypeface(font);
    }

    public void clear(View v){
        EditText name=(EditText) findViewById(R.id.edt_entername);
        name.setText("");
    }

    public void reset(View v){
        TextView txv=(TextView) findViewById(R.id.txv_msg);
        txv.setTextSize(34);
        txv.setTextColor(this.getResources().getColor(R.color.black));
        txv.setText("Hello World!");
        Typeface font=Typeface.create(Typeface.SANS_SERIF,Typeface.NORMAL);
        txv.setTypeface(font);
    }
}