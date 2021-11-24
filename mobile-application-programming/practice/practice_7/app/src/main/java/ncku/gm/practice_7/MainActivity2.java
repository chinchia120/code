package ncku.gm.practice_7;

import androidx.appcompat.app.AppCompatActivity;

import android.content.Intent;
import android.net.Uri;
import android.os.Bundle;
import android.view.View;
import android.widget.Button;
import android.widget.EditText;

public class MainActivity2 extends AppCompatActivity
        implements View.OnClickListener {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main2);

        ((Button)findViewById(R.id.btn_search)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_ok)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_cancel)).setOnClickListener(this);
    }


    @Override
    public void onClick(View view) {
        if(view.getId()==R.id.btn_search){
            startActivity(new Intent(Intent.ACTION_VIEW, Uri.parse("https://www.google.com")));
        }else if(view.getId()==R.id.btn_ok){
            setResult(-1,new Intent().putExtra("店名",((EditText)findViewById(R.id.edt_enter_res)).getText().toString()));
            finish();
        }else if(view.getId()==R.id.btn_cancel){
            finish();
        }
    }
}