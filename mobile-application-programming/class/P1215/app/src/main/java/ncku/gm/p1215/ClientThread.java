package ncku.gm.p1215;

import android.os.Handler;
import android.os.Looper;
import android.os.Message;
import android.util.Log;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.Socket;

public class ClientThread implements Runnable
{
    private Socket mSocket;
    private BufferedReader mBufferedReader = null;
    private OutputStream mOutputStream = null;
    private String host = null;
    private int port = 0;

    private Handler mInputHandler;
    public Handler mOutputHandler;

    public ClientThread(Handler handler, String host, int port)
    {
        this.mInputHandler = handler;
        this.host = host;
        this.port = port;
    }

    @Override
    public void run()
    {
        try
        {
            mSocket = new Socket(host, port);
            Log.d("luhc","connect success "+host+" ("+port+")");
            mBufferedReader = new BufferedReader(new InputStreamReader(mSocket.getInputStream()));
            mOutputStream = mSocket.getOutputStream();

            new Thread()
            {
                @Override
                public void run()
                {
                    super.run();
                    try
                    {
                        String content = null;
                        while ((content = mBufferedReader.readLine()) != null)
                        {
                            Log.d("luhc",content);
                            Message msg = new Message();
                            msg.what = 0;
                            msg.obj = content;
                            mInputHandler.sendMessage(msg);
                        }
                    }
                    catch (IOException e)
                    {
                        e.printStackTrace();
                    }
                }
            }.start();

            Looper.prepare();
            mOutputHandler = new Handler()
            {
                @Override
                public void handleMessage(Message msg)
                {
                    if (msg.what == 1)
                    {
                        try
                        {
                            mOutputStream.write((msg.obj.toString() + "\r\n").getBytes("UTF-8"));
                        }
                        catch (IOException e)
                        {
                            e.printStackTrace();
                        }
                    }
                }
            };
            Looper.loop();
        }
        catch (IOException e)
        {
            e.printStackTrace();
        }
    }
}
