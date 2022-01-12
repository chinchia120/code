package ncku.gm.final_project;

import androidx.appcompat.app.AppCompatActivity;

import android.os.Bundle;
import android.view.View;
import android.widget.Button;
import android.widget.ImageView;
import android.widget.TextView;

public class MainActivity_user_data extends AppCompatActivity {

    UserInformation userInformation = new UserInformation();

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main_user_data);


        ((TextView)findViewById(R.id.txv_show_user_name)).setText("姓名 : "+userInformation.getUser_name());
        ((TextView)findViewById(R.id.txv_show_user_email)).setText("電子郵件 : "+userInformation.getUser_email());
        ((TextView)findViewById(R.id.txv_show_user_password)).setText("密碼 : "+userInformation.getUser_password());
        ((TextView)findViewById(R.id.txv_show_user_phone)).setText("行動電話 : "+userInformation.getUser_phone());

        ((ImageView)findViewById(R.id.imv_back_user_data)).setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                finish();
            }
        });

    }
}