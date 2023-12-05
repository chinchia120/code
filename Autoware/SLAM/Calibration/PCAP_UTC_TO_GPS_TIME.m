%==========================================================================
% Function: Convert UTC time to GPS Week Seconds
% Developed & Coded by: Surachet Srinara, PhD candidate
% Date: July 30, 2021 - 08:49PM
% Copyright@Chet2021
%==========================================================================
function GPS_TIME = PCAP_UTC_TO_GPS_TIME(UTCTimeIn, CFG)

t1_UTCstamp = UTCTimeIn;
t1_UTCadjust = t1_UTCstamp + ((CFG.HOUR-CFG.TimeZone_NOW)*CFG.SEC_IN_ONE_HOUR)*1000000 + CFG.UTC_DELAY;
t1_GPSstamp = (t1_UTCadjust/1000000) + CFG.DOW*CFG.SEC_IN_ONE_DAY + CFG.GPS_LEAP_SEC;

GPS_TIME = t1_GPSstamp;

end