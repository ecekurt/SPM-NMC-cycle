Tamb=[298 310];
temp_fan_on=308;
temp_fan_off=303;

fan_onoff_flag=0;

for d=1:length(Tamb)
temp(1)=Tamb(d);
for k=1:length(power)
    
if  temp(k)>=temp_fan_on
    fan_onoff_flag=1;
    Qohmic =  cur.^2*R_tot_n;
    Qremv=p.h*p.A*(T(k)-p.T_amb);   
%     temp_diff=(heat_gen-heat_rem)/(m_mod*cp_mod); %temperature change in module per each time step
elseif temp(k)<temp_fan_on && temp(k)>=temp_fan_off && fan_onoff_flag_init==1
    fan_onoff_flag=1;
    [heat_gen]=heatgenerated(I,res);
    [heat_rem]=air_cooling(temp_init);
%     temp_diff=(heat_gen-heat_rem)/(m_mod*cp_mod); %temperature change in module per each time step
else
    fan_onoff_flag=0;
    [heat_gen]=heatgenerated_ECM3(I,v_oc,vT);
  heat_rem=0;
   
   temp_diff=(heat_gen-heat_rem)/(m_mod*cp_mod*N_mod); %temperature change in module per each time step
   temp(k+1)=temp(k)+temp_diff;
end
end
end