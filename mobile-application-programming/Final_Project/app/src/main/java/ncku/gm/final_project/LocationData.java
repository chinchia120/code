package ncku.gm.final_project;

public class LocationData {
    private static String dis;
    private static double lat_start,lon_start,lat_end,lon_end;

    public static void setDis(String dis) {
        LocationData.dis = dis;
    }

    public static void setLat_start(double lat_start) {
        LocationData.lat_start = lat_start;
    }

    public static void setLon_start(double lon_start) {
        LocationData.lon_start = lon_start;
    }

    public static void setLon_end(double lon_end) {
        LocationData.lon_end = lon_end;
    }

    public static void setLat_end(double lat_end) {
        LocationData.lat_end = lat_end;
    }

    public static String getDis() {
        return dis;
    }

    public static double getLat_start() {
        return lat_start;
    }

    public static double getLon_start() {
        return lon_start;
    }

    public static double getLat_end() {
        return lat_end;
    }

    public static double getLon_end() {
        return lon_end;
    }
}
