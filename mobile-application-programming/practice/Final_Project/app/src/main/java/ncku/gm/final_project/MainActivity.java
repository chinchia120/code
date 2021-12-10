package ncku.gm.final_project;

import androidx.annotation.Nullable;
import androidx.appcompat.app.AppCompatActivity;

import android.content.Intent;
import android.os.Bundle;
import android.view.View;
import android.widget.Button;
import android.widget.EditText;

import com.google.android.material.snackbar.Snackbar;

import java.util.ArrayList;

public class MainActivity extends AppCompatActivity implements View.OnClickListener {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        ((Button)findViewById(R.id.btn_login)).setOnClickListener(this);
        ((Button)findViewById(R.id.btn_register)).setOnClickListener(this);
    }

    @Override
    public void onClick(View view) {
        int flag=0;
        if(view.getId()==R.id.btn_login){
            for(int i=0;i<str_name.size();i++){
                if(((EditText)findViewById(R.id.edt_enter_email)).getText().toString()==str_name.get(i) || ((EditText)findViewById(R.id.edt_enter_password)).getText().toString()==str_password.get(i)){
                    flag=1;
                    //startActivity(new Intent(this,MainActivity_home.class));
                }
            }
            if(flag==0){
                Snackbar.make(findViewById(R.id.root_login),"登入失敗",Snackbar.LENGTH_LONG).show();
            }else{
                startActivity(new Intent(this,MainActivity_home.class));
            }
        }else if(view.getId()==R.id.btn_register){
            startActivityForResult(new Intent(this,MainActivity_register.class),000);
        }
    }

    ArrayList<String> str_name = new ArrayList<>();
    ArrayList<String> str_email = new ArrayList<>();
    ArrayList<String> str_password = new ArrayList<>();
    ArrayList<String> str_phone_number = new ArrayList<>();

    @Override
    protected void onActivityResult(int requestCode, int resultCode, @Nullable Intent data) {
        super.onActivityResult(requestCode, resultCode, data);

        if(requestCode==000 && resultCode==-1){
            Snackbar.make(findViewById(R.id.root_login),"註冊成功",Snackbar.LENGTH_LONG).show();
            str_name.add(str_name.size(),data.getStringExtra("name"));
            str_email.add(str_email.size(),data.getStringExtra("email"));
            str_password.add(str_password.size(),data.getStringExtra("password"));
            str_phone_number.add(str_phone_number.size(),data.getStringExtra("phone_number"));
        }
    }
}