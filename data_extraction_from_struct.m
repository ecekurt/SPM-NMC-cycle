
oneyear=years(1).days(1:360);
%%
Closs=[];

% for i=1:size(oneyear,2)
% Closs=[Closs;oneyear(i).DChr.CLoss;oneyear(i).Chr.CLoss];
% end
for i=1:size(years,2)
for k=1:size(years(i).days,2)
    Closs=[Closs;years(i).days(k).DChr.CLoss;years(i).days(k).Chr.CLoss];
end
end

%%

Tcell=[];
% for i=1:size(oneyear,2)
% Tcell=[Tcell; oneyear(i).DChr.Tcell; oneyear(i).Chr.Tcell];
% end
for i=1:size(years,2)
for k=1:size(years(i).days,2)
    Tcell=[Tcell;years(i).days(k).DChr.Tcell;years(i).days(k).Chr.Tcell];
end
end
%%

cur=[];
% for i=1:size(oneyear,2)
% Tcell=[Tcell; oneyear(i).DChr.Tcell; oneyear(i).Chr.Tcell];
% end
for i=1:size(years,2)
for k=1:size(years(i).days,2)
     cur=[cur;years(i).days(k).DChr.cur';years(i).days(k).Chr.cur'];
end
end

%%
temp=Tcell;
cap_loss_init=0;
I=cur;

t=0:length(temp);
t_h=t/3600;
% del_Ih=zeros(1,length(temp));
% d_cyc_loss=zeros(1,length(temp));
cyc_loss=zeros(1,length(temp));
cyc_loss(1)=cap_loss_init;

for i=2:length(temp)
    
%     
%     Crate(i)=abs(I(i-1))/2.3;
%     if Crate(i)>Crate_index(end)
%         %         A(i)=A_map(end);
%         Crate(i)=Crate_index(end);
%     elseif Crate(i)<Crate_index(1)
%         %         A(i)=A_map(1);
%         Crate(i)=Crate_index(1);
%         
%     end
%     A=2.95e+05;
%     B=7354;
%     n=1.285;

    A=0.008329;
    B=713.6;
    n=0.62;

    Ih_in=(cyc_loss(i-1)/(A*exp(-B/((temp(i))))))^(1/n);
    %     del_Ih(i)=abs(I(i-1))*(t_h(i)-t_h(i-1)); %Ah processed per cell
    delt=t_h(i)-t_h(i-1);
    if I(i)*I(i-1)<0
        
        a=abs(I(i-1))*delt/(abs(I(i))+abs(I(i-1)));
        del_Ih(i)=abs(I(i-1))*a/2+abs(I(i))*(delt-a)/2;
    else
        del_Ih(i)=(abs(I(i))+abs(I(i-1)))*delt/2;
    end
    d_cyc_loss(i)=A*exp(-B/((temp(i))))*((Ih_in+del_Ih(i))^n-Ih_in^n);
    cyc_loss(i)=cyc_loss(i-1)+d_cyc_loss(i);
    
end
cap_loss=cyc_loss(end);