package ncku.gm.final_project;

import androidx.appcompat.app.AppCompatActivity;

import android.content.Intent;
import android.os.Bundle;

import android.os.Handler;
import android.os.Message;

import android.view.View;
import android.widget.Button;
import android.widget.ImageView;
import android.widget.EditText;
import android.widget.TextView;

import com.google.android.material.snackbar.Snackbar;

import java.util.Random;

public class MainActivity_chat extends AppCompatActivity implements View.OnClickListener {

    ClientThread mClientThread;
    private Handler mInputHandler;
    private String host = "140.116.47.94";
    private int port = 7070;
    UserInformation userInformation = new UserInformation();


    //??????
    TextView tvmessage;
    ImageView imgsend;
    EditText edtxmessage;
    static String chat;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main_chat);

        ((ImageView)findViewById(R.id.imgsend)).setOnClickListener(this);
        ((ImageView)findViewById(R.id.imv_back_chat)).setOnClickListener(this);

        edtxmessage = findViewById(R.id.message);
        imgsend = findViewById(R.id.imgsend);
        tvmessage = findViewById(R.id.chatbox);

        ((TextView)findViewById(R.id.chatbox)).setText(chat);

        mInputHandler = new Handler()
        {
            @Override
            public void handleMessage(Message msg)
            {
                if (msg.what == 0)
                {
                    String[] tmp = msg.obj.toString().split(",");
                    //???????
                    if(tmp[0].matches("message")){

                        //Snackbar.make(findViewById(R.id.root_chat),"Get Message",Snackbar.LENGTH_LONG).show();
                        //???????????????X?????????(???????,??????,?????????)
                        /*
                        if( (tmp[1].matches(userInformation.getUser_name()) || tmp[3].matches(userInformation.getUser_name()) || (tmp[1].matches(tmp[3])))){
                            tvmessage.append(String.format("%s to %s : %s\n",tmp[1],tmp[3],tmp[2]));
                            chat = tvmessage.getText().toString();
                        }

                         */
                        tvmessage.append(String.format("%s : %s\n",tmp[1],tmp[2]));
                        chat = tvmessage.getText().toString();
                    }
                }
            }
        };

        mClientThread = new ClientThread(mInputHandler, host, port);
        new Thread(mClientThread).start();
    }

    @Override
    public void onClick(View view) {
        if (view.getId() == R.id.imgsend){
            if(!edtxmessage.getText().toString().matches("")){
                Message msg = new Message();
                msg.what = 1;
                msg.obj = "message" + "," +
                        userInformation.getUser_name()+","+
                        edtxmessage.getText().toString()+","+
                        getIntent().getBundleExtra("Bundle").getString("name");;
                mClientThread.mOutputHandler.sendMessage(msg);
                edtxmessage.setText("");
            }
        }
        else if (view.getId() == R.id.imv_back_chat){
            finish();
        }
    }
}