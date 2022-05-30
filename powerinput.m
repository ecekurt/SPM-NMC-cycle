function [current,C_rate] = powerinput(data,V0,t,p)

power= interp1(data.power_index,data.power,t); 

current=power/V0;

C_rate=current./p.Ah;

end