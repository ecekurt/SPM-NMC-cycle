function [Qohmic,Qremv]=coolingcontrol(T,cur,R_tot_n,p)
temp_fan_on=308;
temp_fan_off=303;

fan_onoff_flag=0;

if  T>=temp_fan_on
    fan_onoff_flag=1;
    Qohmic =  cur.^2*R_tot_n;
    Qremv= p.h*p.A*(T-p.T_amb);   

elseif T<temp_fan_on && T>=temp_fan_off && fan_onoff_flag==1
    fan_onoff_flag=1;
    Qohmic =  cur.^2*R_tot_n;
    Qremv=p.h*p.A*(T-p.T_amb); 

else
    fan_onoff_flag=0;
    Qohmic =  cur.^2*R_tot_n;
    Qremv=0;
   

end
