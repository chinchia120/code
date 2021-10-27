package ncku.gm.p1027;

import androidx.appcompat.app.AppCompatActivity;

import android.os.Bundle;
import android.view.View;
import android.widget.AdapterView;
import android.widget.ArrayAdapter;
import android.widget.Button;
import android.widget.ListView;
import android.widget.Spinner;
import android.widget.TextView;

import java.util.ArrayList;

public class MainActivity extends AppCompatActivity
        implements View.OnClickListener , AdapterView.OnItemSelectedListener , AdapterView.OnItemClickListener {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        ((Button)findViewById(R.id.btn_book)).setOnClickListener(this);
        ((Spinner)findViewById(R.id.spn_theater)).setOnItemSelectedListener(this);
        ((ListView)findViewById(R.id.lsv_food)).setOnItemClickListener(this);

        for(int i=0;i<4;i++){
            order[i]=false;
        }

        String[] mt = {"台北影城","台中影城","台南影城","高雄影城"};
        ArrayAdapter<String> ad = new ArrayAdapter<>(this, android.R.layout.simple_spinner_dropdown_item,mt);
        ((Spinner)findViewById(R.id.spn_theater)).setAdapter(ad);
        //ad.setDropDownViewResource(android.R.layout.simple_spinner_dropdown_item);

        ArrayList<String> fd2 = new ArrayList<>();
        fd2.add("漢堡");

        String[] fd = {"漢堡","薯條","可樂","濃湯","炸雞"};
        ArrayAdapter<String> ad2 = new ArrayAdapter<>(this, android.R.layout.simple_list_item_1,fd);
        ((ListView)findViewById(R.id.lsv_food)).setAdapter(ad2);
    }

    @Override
    public void onClick(View view) {
        show();
    }

    public void show(){
        Spinner spn = findViewById(R.id.spn_theater);
        /*
        int idx = spn.getSelectedItemPosition();
        if(idx==0){
            ((TextView)findViewById(R.id.txv_msg)).setText("訂購台北影城的票");
        }else if(idx==1){
            ((TextView)findViewById(R.id.txv_msg)).setText("訂購台中影城的票");
        }else{
            ((TextView)findViewById(R.id.txv_msg)).setText("訂購台南影城的票");
        }
        */
        String str = spn.getSelectedItem().toString();
        ((TextView)findViewById(R.id.txv_msg)).setText("訂購" + str + "的票");
    }

    @Override
    public void onItemSelected(AdapterView<?> adapterView, View view, int position, long id) {
        /*
        if(position==0){
            ((TextView)findViewById(R.id.txv_msg)).setText("訂購台北影城的票");
        }else if(position==1){
            ((TextView)findViewById(R.id.txv_msg)).setText("訂購台中影城的票");
        }else{
            ((TextView)findViewById(R.id.txv_msg)).setText("訂購台南影城的票");
        }

        String str = ((TextView)view).getText().toString();
        ((TextView)findViewById(R.id.txv_msg)).setText("訂購" + str + "的票");

        show();
        */
    }

    @Override
    public void onNothingSelected(AdapterView<?> adapterView) {

    }
    boolean[] order = new boolean[4];

    @Override
    public void onItemClick(AdapterView<?> adapterView, View view, int position, long id) {
        /*
        Spinner spn = findViewById(R.id.spn_theater);
        String str_theater = spn.getSelectedItem().toString();


        if(order[position]==false) {
            order[position]=true;
        }else{
            order[position]=false;
        }

        String str="你點了";
        for(int i=0;i<4;i++){
            if(order[i]==true){
                if(i==0){
                    str+="漢堡";
                }else if(i==1){
                    str+="薯條";
                }else if(i==2){
                    str+="可樂";
                }else {
                    str+="玉米濃湯";
                }
            }
        }

        ((TextView)findViewById(R.id.txv_msg)).setText("你點了 : "  + str );
        */
        /*
        if(position==0){
            ((TextView)findViewById(R.id.txv_msg)).setText("你點了 : 漢堡");
        }else if(position==1){
            ((TextView)findViewById(R.id.txv_msg)).setText("你點了 : 薯條");
        }else if(position==2){
            ((TextView)findViewById(R.id.txv_msg)).setText("你點了 : 可樂");
        }else{
            ((TextView)findViewById(R.id.txv_msg)).setText("你點了 : 玉米濃湯");
        }

        String str_food = ((TextView)view).getText().toString();
        ((TextView)findViewById(R.id.txv_msg)).setText("你點了 : "  + str_food );
        */
        /*
        ArrayList <String> select = new ArrayList<>();
        String str =((TextView)view).getText().toString();
        if(select.contains(str)){
            select.remove(str);
        }else{
            select.add(str);
        }
        String str2 ="你點了 : ";
        for(int i=0;i<select.size();i++){
            str2 +=  select.get(i) + " ";
        }
        for(String s : select){
            str2 += s + " ";
        }
        ((TextView)findViewById(R.id.txv_msg)).setText(str2);
        */

    }

}