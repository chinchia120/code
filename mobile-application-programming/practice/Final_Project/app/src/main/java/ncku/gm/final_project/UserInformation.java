package ncku.gm.final_project;

public class UserInformation {
    private static String user_name,user_email,user_password,user_phone;

    public static void setUser_name(String user_name) {
        UserInformation.user_name = user_name;
    }

    public static void setUser_email(String user_email) {
        UserInformation.user_email = user_email;
    }

    public static void setUser_password(String user_password) {
        UserInformation.user_password = user_password;
    }

    public static void setUser_phone(String user_phone) {
        UserInformation.user_phone = user_phone;
    }

    public static String getUser_name() {
        return user_name;
    }

    public static String getUser_phone() {
        return user_phone;
    }

    public static String getUser_email() {
        return user_email;
    }

    public static String getUser_password() {
        return user_password;
    }
}
