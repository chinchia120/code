%% =============== Function of PCAP_UTC_TO_GPS_TIME =============== %%
function GPS_TIME = PCAP_UTC_TO_GPS_TIME(UTCTimeIn, CFG)
    t1_UTCstamp = UTCTimeIn;
    t1_UTCadjust = t1_UTCstamp + ((CFG.HOUR-CFG.TimeZone_NOW)*CFG.SEC_IN_ONE_HOUR)*1000000 + CFG.UTC_DELAY;
    t1_GPSstamp = (t1_UTCadjust/1000000) + CFG.DOW*CFG.SEC_IN_ONE_DAY + CFG.GPS_LEAP_SEC;

    GPS_TIME = t1_GPSstamp;
end