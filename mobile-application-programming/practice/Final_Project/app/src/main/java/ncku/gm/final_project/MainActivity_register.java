package ncku.gm.final_project;

import androidx.appcompat.app.AppCompatActivity;

import android.content.Intent;
import android.os.Bundle;
import android.view.View;
import android.widget.Button;
import android.widget.EditText;

import com.google.android.material.snackbar.Snackbar;

public class MainActivity_register extends AppCompatActivity implements View.OnClickListener {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main_register);

        ((Button)findViewById(R.id.btn_ok_register)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_cancel_register)).setOnClickListener(this);
    }

    @Override
    public void onClick(View view) {
        if(view.getId()==R.id.btn_ok_register){
            if(((EditText)findViewById(R.id.edt_enter_name)).getText().toString().matches("") || ((EditText)findViewById(R.id.edt_enter_email_register)).getText().toString().matches("") || ((EditText)findViewById(R.id.edt_enter_password_register)).getText().toString().matches("") || ((EditText)findViewById(R.id.edt_enter_phone)).getText().toString().matches("")){
                Snackbar.make(findViewById(R.id.root_register),"資訊不完整",Snackbar.LENGTH_LONG).show();
            }else{
                Intent it = new Intent();
                it.putExtra("Name",((EditText)findViewById(R.id.edt_enter_name)).getText().toString());
                it.putExtra("Email",((EditText)findViewById(R.id.edt_enter_email_register)).getText().toString());
                it.putExtra("Password",((EditText)findViewById(R.id.edt_enter_password_register)).getText().toString());
                it.putExtra("Phone",((EditText)findViewById(R.id.edt_enter_phone)).getText().toString());
                setResult(RESULT_OK,it);
                finish();
            }
        }else if(view.getId()==R.id.btn_cancel_register){
            finish();
        }
    }
}