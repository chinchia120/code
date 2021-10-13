package ncku.gm.p1013;

import androidx.appcompat.app.AppCompatActivity;

import android.content.Context;
import android.os.Bundle;
import android.os.Vibrator;
import android.view.MotionEvent;
import android.view.View;
import android.widget.Button;
import android.widget.TextView;

public class MainActivity extends AppCompatActivity
        implements View.OnClickListener , View.OnLongClickListener , View.OnTouchListener {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        Button btn = findViewById(R.id.btn_plus1);
        btn.setOnClickListener(this);
        btn.setOnLongClickListener(this);

        TextView txv = findViewById(R.id.txv_msg);
        txv.setOnLongClickListener(this);

        TextView txv2 = findViewById(R.id.txv_vbt);
        txv2.setOnTouchListener(this);
    }

    int cnt=0;
    Vibrator vb;
    @Override
    public void onClick(View view) {
        TextView txv = findViewById(R.id.txv_msg);
        txv.setText(++cnt + "");
    }

    @Override
    public boolean onLongClick(View view) {
        TextView txv = findViewById(R.id.txv_msg);
        if(view.getId()==R.id.btn_plus1){
            cnt+=2;
            txv.setText(cnt + "");
        }else{
            cnt=0;
            txv.setText(cnt + "");
        }
        return true;
    }

    @Override
    public boolean onTouch(View view, MotionEvent motionEvent) {
        vb = (Vibrator) getSystemService(Context.VIBRATOR_SERVICE);
        if(motionEvent.getAction()==MotionEvent.ACTION_DOWN){
            vb.vibrate(5000);
        }else if(motionEvent.getAction()==MotionEvent.ACTION_UP){
            vb.cancel();
        }
        return true;
    }
}