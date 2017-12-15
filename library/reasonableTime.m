function [timeVal, unit] = reasonableTime(timeVal)
unit = 's';
if(strcmp(unit, 's'))
    if(abs(timeVal) > 60)
        unit = 'min';
        timeVal = timeVal/60;
    end
end
if(strcmp(unit, 'min'))
    if(abs(timeVal) > 60)
        unit = 'hrs';
        timeVal = timeVal/60;
    end
end
if(strcmp(unit, 'hrs'))
    if(abs(timeVal) > 24)
        unit = 'days';
        timeVal = timeVal/24;
    end
end
if(strcmp(unit, 'days'))
    if(abs(timeVal) > 180)
        unit = 'years';
        timeVal = timeVal/365.242199;
    end
end
return