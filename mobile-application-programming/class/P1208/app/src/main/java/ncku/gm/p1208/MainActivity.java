package ncku.gm.p1208;

import androidx.annotation.Nullable;
import androidx.appcompat.app.AppCompatActivity;
import androidx.core.app.ActivityCompat;

import android.Manifest;
import android.app.Activity;
import android.content.ContentValues;
import android.content.Intent;
import android.content.pm.ActivityInfo;
import android.content.pm.PackageManager;
import android.graphics.Bitmap;
import android.graphics.BitmapFactory;
import android.graphics.Matrix;
import android.net.Uri;
import android.os.Bundle;
import android.provider.MediaStore;
import android.view.View;
import android.widget.ImageView;
import android.widget.Toast;

public class MainActivity extends AppCompatActivity {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        if(ActivityCompat.checkSelfPermission(this, Manifest.permission.WRITE_EXTERNAL_STORAGE)!= PackageManager.PERMISSION_GRANTED){
            ActivityCompat.requestPermissions(this,new String[]{Manifest.permission.WRITE_EXTERNAL_STORAGE},200);
        }

        setRequestedOrientation(ActivityInfo.SCREEN_ORIENTATION_NOSENSOR);
        setRequestedOrientation(ActivityInfo.SCREEN_ORIENTATION_PORTRAIT);
    }

    Uri uri;

    public void take_picture(View v){
        uri = getContentResolver().insert(MediaStore.Images.Media.EXTERNAL_CONTENT_URI,new ContentValues());
        Intent it = new Intent(MediaStore.ACTION_IMAGE_CAPTURE);
        it.putExtra(MediaStore.EXTRA_OUTPUT,uri);
        startActivityForResult(it,100);
        //startActivityForResult(new Intent(MediaStore.ACTION_IMAGE_CAPTURE),100);

        it.putExtra(MediaStore.EXTRA_OUTPUT,uri);

    }

    @Override
    protected void onActivityResult(int requestCode, int resultCode, @Nullable Intent data) {
        super.onActivityResult(requestCode, resultCode, data);

        if(requestCode==100 && resultCode== Activity.RESULT_OK){
            //Bundle bdl = data.getExtras();
            //Bitmap bmp = (Bitmap) bdl.get("data");
            //((ImageView)findViewById(R.id.imv_show)).setImageBitmap((Bitmap)data.getExtras().get("data"));
            Bitmap bmp = null;
            try {
                BitmapFactory.Options opt = new BitmapFactory.Options();
                opt.inJustDecodeBounds = true;
                bmp = BitmapFactory.decodeStream(getContentResolver().openInputStream(uri),null,opt);
                opt.inSampleSize = 2;
                opt.inJustDecodeBounds = false;
                bmp = BitmapFactory.decodeStream(getContentResolver().openInputStream(uri),null,opt);

                Matrix m = new Matrix();
                m.postRotate(90);
                bmp = Bitmap.createBitmap(bmp,0,0,opt.outWidth,opt.outHeight,m,true);

                Toast.makeText(this,opt.outWidth+","+opt.outHeight,Toast.LENGTH_LONG).show();
            }catch (Exception e){

            }
            ((ImageView)findViewById(R.id.imv_show)).setImageBitmap(bmp);
        }else if(requestCode==101 && resultCode==Activity.RESULT_OK){
            Uri uri = data.getData();
            Bitmap bmp = null;
            try {
                bmp = BitmapFactory.decodeStream(getContentResolver().openInputStream(uri),null,null);
            }catch (Exception e){

            }
            ((ImageView)findViewById(R.id.imv_show)).setImageBitmap(bmp);
        }
    }

    public void pick_picture(View v){
        Intent it = new Intent(Intent.ACTION_GET_CONTENT);
        it.setType("image/*");
        startActivityForResult(it,101);
    }
}