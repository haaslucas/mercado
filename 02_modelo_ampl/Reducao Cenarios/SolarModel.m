%
%   SOLAR MODEL BASED ON "An Introduction to Environmental Biophysics: 2nd edition - Gaylon S. Campbell, ?John M. Norman (1998)"
%
%   INPUTS: - Latitude in degrees
%           - Airpressure in the Altitude considered in Pascal (Pa)
%
%   OUTPUT: - Monthly Average Irradiation (W/m�)
%           - Daily Average Irradiation (W/m�)
%
function [Monthly,Hourly]=SolarModel(Latitude, airpressure)

    Month=zeros(1,12);       
    p=airpressure; 
    Hourly=0;
    
    for DayOfYear=1:365
            t=0.75; % transmittance (unitless)
            S=1367; %solar constant (w/m^2)

            hours=[7,8,9,10,11,12,13,14,15,16,17];
            hangle=(12.0-hours)*15.0*pi/180.0; % Note pi/180.0 is the factor to convert angles in degrees in radiances.

            declangle=23.45*sin(2.0*pi*(284.0+DayOfYear)/365.0)*pi/180.0;

            cosz=sin(Latitude*pi/180.0)*sin(declangle)+cos(Latitude*pi/180.0)*cos(declangle)*cos(hangle);

            m=p./(101.3*cosz); %optical airmass

            Sb=cosz*S.*(t.^m); % Sb is the beam radiation on a horizontal surface
            Sd=0.3*(1.0-t.^m)*S.*cosz; % Sd is the diffuse radiation on a horizontal surface
            St=Sb+Sd; % St is the total radiation;
            
            Hourly=(Hourly+St);
            Irradiance_Wm2 = mean(St);

            %Jan = 1 - 31 (31)
            %Feb = 32 - 59 (28)
            %Mar = 60 - 90 (31)
            %Abr = 91 - 120 (30)
            %May = 121 - 151 (31)
            %Jun = 152 - 181 (30)
            %Jul = 182 - 212 (31)
            %Aug = 213 - 243 (31)
            %Sep = 244 - 273 (30)
            %Oct = 274 - 304 (31)
            %Nov = 305 - 334 (30)
            %Dec = 335 - 365 (31)

            if DayOfYear <= 31
                Month(1,1) = Month(1,1)+Irradiance_Wm2;
            elseif DayOfYear >=32 && DayOfYear <=59
                Month(1,2) = Month(1,2)+Irradiance_Wm2;
            elseif DayOfYear >=60 && DayOfYear <=90
                Month(1,3) = Month(1,3)+Irradiance_Wm2;
            elseif DayOfYear >=91 && DayOfYear <=120
                Month(1,4) = Month(1,4)+Irradiance_Wm2;
            elseif DayOfYear >=121 && DayOfYear <=151
                Month(1,5) = Month(1,5)+Irradiance_Wm2;
            elseif DayOfYear >=152 && DayOfYear <=181
                Month(1,6) = Month(1,6)+Irradiance_Wm2;
            elseif DayOfYear >=182 && DayOfYear <=212
                Month(1,7) = Month(1,7)+Irradiance_Wm2;
            elseif DayOfYear >=213 && DayOfYear <=243
                Month(1,8) = Month(1,8)+Irradiance_Wm2;
            elseif DayOfYear >=244 && DayOfYear <=273
                Month(1,9) = Month(1,9)+Irradiance_Wm2;
            elseif DayOfYear >=274 && DayOfYear <=304
                Month(1,10) = Month(1,10)+Irradiance_Wm2;
            elseif DayOfYear >=305 && DayOfYear <=334
                Month(1,11) = Month(1,11)+Irradiance_Wm2;
            elseif DayOfYear >=335 && DayOfYear <=365
                Month(1,12) = Month(1,12)+Irradiance_Wm2;
            end
    end
    Month(1,1)=Month(1,1)/31;
    Month(1,2)=Month(1,2)/28;
    Month(1,3)=Month(1,3)/31;
    Month(1,4)=Month(1,4)/30;
    Month(1,5)=Month(1,5)/31;
    Month(1,6)=Month(1,6)/30;
    Month(1,7)=Month(1,7)/31;
    Month(1,8)=Month(1,8)/31;
    Month(1,9)=Month(1,9)/30;
    Month(1,10)=Month(1,10)/31;
    Month(1,11)=Month(1,11)/30;
    Month(1,12)=Month(1,12)/31;
    Monthly=Month;
    Hourly=Hourly/365;
end
