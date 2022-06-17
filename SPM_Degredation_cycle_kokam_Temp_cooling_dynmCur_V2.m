%% KOKAM NMC battery SPM model with SEI effect
% SLIDE Kokam NMC parameters
% Every cycle starts with a different ambient temperature.
% Cooling control was added. Temperature is controled by adjusting the convection coefficient h[W/(m2 K)]. 


clear all
close all
clc 

run Kokam_NMC_parameters
load('KokamOCVNMC.mat'); %half cell pot.
load('KokamOCVC.mat');  %half cell pt.
load('OCVcell.mat');  %overall cell ocv vs Ah
load('KokamNMC.mat'); 
load('KokamC.mat');   
load('DailyTemp_PHOENIX.mat');
load('cur2.mat'); % dynamic input current
load('cur.mat'); % dynamic input current
load('Pbatt.mat'); 
global data 
global p
global KokamOCVNMC KokamNMC 
global KokamOCVC KokamC
global OCVcell



% Charge Input

p.C_rate= 0.556;
p.cycle_number=2;
p.Nc= p.cycle_number*2;

tend_chr= 3600*p.Nc/p.C_rate;
data.chrcur=@(t) 0*(t==0) - 2.7*p.C_rate*(t<=1) + 2.7*p.C_rate*(1<=tend_chr); % charge current - constant charge
data.chrtime=1:tend_chr;


% Discharge Input
data.power=data.power;
data.power_index=data.power_index;
tend= length(data.power_index);
V0=4.2;

%% Finite difference for spherical particle and electrolyte
p.Np=50;
p.Nn=50;
p.delta_p =  p.R_p/(p.Np);
p.delta_n =  p.R_n/(p.Nn);  

%% Initial conditions 

% Initial concentrations of solid particles

theta_p_min=p.theta_p_min;
theta_n_max=p.theta_n_max;
Up0=p.c_s_p_max*theta_p_min*ones(p.Np-1,1);             
Un0=p.c_s_n_max*theta_n_max*ones(p.Nn-1,1);            

% Temperature
T10 = 298.15; %Core Temp.
T20 = 298.15; %Surface Temp.
DailyT= [298.15 320.15]; 
% DailyT= DailyTemp; % Real Temp. data for 366 days

% SEI
Qs0=p.eps_s_n*p.Faraday*p.Area_n*p.L_n*p.c_s_n_max*p.theta_n_max/1.0844; %scaling to 2.7
sei0=p.L_sei;

tic  
 
% set up starting equation    
tspan=0:1:tend;
t0=tspan(1);
x0 = [Un0; Up0; T10; T20; DailyT(1);Qs0; sei0]';
func=@ode_SPMT_discharge;
options=odeset('Events',@Efcn); 
year=1;

for i=1:year

    for j=1:length(DailyT)
    
    % Update the temperature
    T0=DailyT(j);
    x0(:,end-2)=T0;

    % Run integration until event function stops it
    options=odeset('Events',@Efcn); 
    [tDChr,xDChr] = ode23s(@(t,x) ode_SPMT_discharge(t,x,V0),[1:tend-1],x0,options); 
    
    % Check the concentrations at the surface of the particles     
    if (xDChr(end,49) <= p.c_s_n_max*p.theta_n_min || xDChr(end,49) > p.c_s_n_max*p.theta_n_min  || xDChr(end,98) >= p.c_s_p_max*p.theta_p_max)
    
        options=odeset('Events',@Efcn1);    
        [tChr,xChr] = ode23s(@ode_SPMT_charge,[1:tend_chr],xDChr(end, :),options); 
        x0=xChr(end,:);
        flag=1
    
    end
     
     % Store the data
            for k=1:length(tDChr)
                [~,theta_p(k),theta_n(k),V_spm(k),V_ocv(k), eta_n(k), eta_p(k), ...
                     Uref_n(k), Qohmic(k),Qremv(k),R_tot_n(k),cur(k)]...
                    =ode_SPMT_discharge(tDChr(k),xDChr(k,:)',V0);
                SOCp(k)=( theta_p(k)- p.theta_p_max )/( p.theta_p_min -p.theta_p_max);
                SOCn(k)=( theta_n(k)- p.theta_n_min )/( p.theta_n_max -p.theta_n_min);
        
            end
            
            years(i).days(j).DChr.Voltage=V_spm;
            years(i).days(j).DChr.NLicon_n=theta_n;
            years(i).days(j).DChr.NLicon_p=theta_p;
            years(i).days(j).DChr.SOCn=SOCn;
            years(i).days(j).DChr.SOCp=SOCp;
            years(i).days(j).DChr.cur=cur;
            years(i).days(j).DChr.Tcell=xDChr(:,end-2);
            years(i).days(j).DChr.cn=xDChr(:,(p.Nn-1));
            years(i).days(j).DChr.cp=xDChr(:,2*(p.Nn-1));
            years(i).days(j).DChr.CLoss=xDChr(:,end-1)./3600;
            years(i).days(j).DChr.SEIgrowth=xDChr(:,end);
    
            clear VChr_spm thetaChr_p thetaChr_n SOCn SOCp
               for k=1:length(tChr)
            [~,thetaChr_p(k),thetaChr_n(k),VChr_spm(k),V_ocv(k), eta_n(k), eta_p(k), ...
                     Uref_n(k), Qohmic(k),Qremv(k),R_tot_n(k),cur(k)]...
                =ode_SPMT_charge(tChr(k),xChr(k,:)');
            SOCp(k)=( thetaChr_p(k)- p.theta_p_max )/( p.theta_p_min -p.theta_p_max);
            SOCn(k)=( thetaChr_n(k)- p.theta_n_min )/( p.theta_n_max -p.theta_n_min);
               end  
    
            years(i).days(j).Chr.Voltage=VChr_spm;
            years(i).days(j).Chr.NLicon_n=thetaChr_n;
            years(i).days(j).Chr.NLicon_p=thetaChr_p;
            years(i).days(j).Chr.SOCn=SOCn;
            years(i).days(j).Chr.SOCp=SOCp;
            years(i).days(j).Chr.cur=cur;
            years(i).days(j).Chr.Tcell=xChr(:,end-2);
            years(i).days(j).Chr.cn=xChr(:,(p.Nn-1));
            years(i).days(j).Chr.cp=xChr(:,2*(p.Nn-1));
            years(i).days(j).Chr.CLoss=xChr(:,end-1)./3600;
            years(i).days(j).Chr.SEIgrowth=xChr(:,end);
                  
      
    end   
    
   % Update the initial states 
    x0=xChr(end,:);
    V0=max(VChr_spm);
end

toc




function [xdot,varargout]=ode_SPMT_discharge(t,x,V0)

global data
global p
global KokamOCVNMC KokamNMC 
global KokamOCVC KokamC
global OCVcell

U_n = x(1:(p.Nn-1));
U_p = x(p.Nn : 2*(p.Nn-1));
T1 = real(x(end-4));
T2 = real(x(end-3));
T  = real(x(end-2));
Q_s= real(x(end-1));
sei = real(x(end));


[cur,p.C_rate] = powerinput(data,V0,t,p);


TEMP=T;

%% Solid phase dynamics

% Molar flux for solid phase [mol m-2 s-1]
J_p=-(cur./p.Area_p)./(p.Faraday*p.a_p*p.L_p); 
J_n=(cur./p.Area_n)./(p.Faraday*p.a_n*p.L_n);  

% Solid phase diffusivity temperature dependence [m s-1]
p.Ds_n = p.Ds_n0 * exp(p.E.Dsn/p.R*(1/p.T_ref - 1/TEMP)); 
p.Ds_p = p.Ds_p0 * exp(p.E.Dsp/p.R*(1/p.T_ref - 1/TEMP)); 

% Matrices for solid-phase Li concentration
 [A_p,A_n,B_n,B_p]= matrixs(p);

% Surface concentrations

c_ss_p= U_p(end);
c_ss_n= U_n(end);
   

%% Calculation of potential of the cell

% li-fraction  of the electrodes  
 [theta_p,theta_n]=getsoc(c_ss_p,c_ss_n,p);
 
% OCV of the cell [v]

[Uref_p, Uref_n]=refpotantial (theta_p, theta_n,KokamOCVNMC, KokamNMC, KokamOCVC, KokamC, OCVcell);
V_ocv = Uref_p-Uref_n  ;

% Kinetic reaction rate, adjusted for Arrhenius temperature dependence [(A/m^2)*(mol^3/mol)^(1+alpha)]

 p.k_n = p.k_n0 * exp(p.E.kn/p.R*(1/p.T_ref - 1/TEMP));   
 p.k_p = p.k_p0 * exp(p.E.kp/p.R*(1/p.T_ref - 1/TEMP)); 

% Exchange current density [I m-2]

 i_0n = p.k_n * p.Faraday * sqrt(((p.c_s_n_max - c_ss_n) .* c_ss_n .* p.ce));
 i_0p = p.k_p * p.Faraday * sqrt(((p.c_s_p_max - c_ss_p) .* c_ss_p .* p.ce));

 xn= 0.5 * (cur/p.Area_n)/(p.a_n*p.L_n)/i_0n;
 xp=-0.5 * (cur/p.Area_p)/(p.a_p*p.L_p)/i_0p;

% Overpotentials [V]
 RTaF=(p.R*TEMP)./(p.alph*p.Faraday);
 
 eta_n  = RTaF .* asinh(cur ./ (2.*p.a_n.*p.Area_n.*p.L_n.*i_0n));
 eta_p  = RTaF .* asinh(-cur ./ (2.*p.a_p.*p.Area_p.*p.L_p.*i_0p));

% eta_p=(2*p.R*TEMP)/p.Faraday * log(xp + sqrt(1+xp*xp));
% eta_n=(2*p.R*TEMP)/p.Faraday * log(xn + sqrt(1+xn*xn));

% SPM Voltage [V]
 V_spm= eta_p - eta_n + V_ocv;
 
% Update the voltage for current calculation
 V0=V_spm;

%% Degredation  


Sn=3*p.eps_s_n*p.L_n*p.Area_n/p.R_n;  %m^2                                         
  
% Overpotential of the SEI [V]

eta_sei_n= Uref_n + eta_n - p.Us + (sei*2037.4)*cur;       %2037.4 R                         

% Current density of the SEI [I m-2]

ksei = p.ksei* exp((130000/p.R)*(1/p.T_ref - 1/TEMP)); 
is= real (p.Faraday*ksei*exp(-p.Faraday/(p.R*TEMP) *( eta_sei_n) ));                        

% Growth rate of SEI layer [m]

sei_dot = is/(p.Faraday*p.rhos);

% Total resistance (film + growing SEI layer)

R_tot_n = p.Rsei_n*sei + p.Rcell;

% cyclable Li, evolution of the charge 

Q_dot= -Sn*is;                                                                 

% Total molar flux [mol m-2 s-1]

js_n = (J_n) + (is/p.Faraday);

% Lithium Loss [A]

LLi=is*p.Area_n*p.L_n*p.a_n;


% Solid particle concentration [mol m-3]

c_p = A_p*U_p + B_p.*J_p;
c_n = A_n*U_n + B_n.*js_n;  

%% Heat generation  

% Heat remove

[Qohmic,Qremv]=coolingcontrol(T,cur,R_tot_n,p);
Qgen= Qohmic + 0 + 0;

% Temperature calculation [K]
 T1_dot= Qgen./p.Cc + (T2-T1)./(p.Rc*p.Cc); %Tc core temp
 T2_dot= (p.T_cool-T2)./(p.Ru*p.Cs) - (T2-T1)./(p.Rc*p.Cs); %Ts surface tem. 

% Lumped Temperature calculation [K]

T_dot= (Qgen  - Qremv)./(p.M*p.Cp);

%% Outputs
xdot = [c_n; c_p; real(T1_dot); real(T2_dot); real(T_dot); real(Q_dot); real(sei_dot)]; 

varargout{1} = theta_p;
varargout{2} = theta_n;
varargout{3} = V_spm;
varargout{4} = V_ocv;
varargout{5} = eta_n;
varargout{6} = LLi;
varargout{7}= Uref_n;
varargout{8}= Qohmic;
varargout{9}= Qremv;
varargout{10}= R_tot_n;
varargout{11}= cur;
end



function [xdot,varargout]=ode_SPMT_charge(t,x)
global data
global p
global KokamOCVNMC KokamNMC 
global KokamOCVC KokamC
global OCVcell


U_n = x(1:(p.Nn-1));
U_p = x(p.Nn : 2*(p.Nn-1));
T1 = real(x(end-4));
T2 = real(x(end-3));
T  = real(x(end-2));
Q_s= real(x(end-1));
sei = real(x(end));


cur=-data.chrcur(t);


TEMP=T;

%% Solid phase dynamics

% Molar flux for solid phase [mol m-2 s-1]
J_p=-(cur./p.Area_p)./(p.Faraday*p.a_p*p.L_p);
J_n=(cur./p.Area_n)./(p.Faraday*p.a_n*p.L_n);

% Solid phase diffusivity temperature dependence [m s-1]
p.Ds_n = p.Ds_n0 * exp(p.E.Dsn/p.R*(1/p.T_ref - 1/TEMP));
p.Ds_p = p.Ds_p0 * exp(p.E.Dsp/p.R*(1/p.T_ref - 1/TEMP)) ;

% Matrices for solid-phase Li concentration
 [A_p,A_n,B_n,B_p]= matrixs(p);

% Calculation of the surface concentration [mol m-3]

c_ss_p= U_p(end);
c_ss_n= U_n(end);
 

%% Calculation of potential of the cell

% li-fraction of the electrodes  
 [theta_p,theta_n]=getsoc(c_ss_p,c_ss_n,p);
  
% OCV of the cell [V]

 [Uref_p, Uref_n]=refpotantial (theta_p, theta_n,KokamOCVNMC, KokamNMC, KokamOCVC, KokamC, OCVcell);
 V_ocv = Uref_p-Uref_n  ;

% Kinetic reaction rate, adjusted for Arrhenius temperature dependence [(A/m^2)*(mol^3/mol)^(1+alpha)]

 p.k_n = p.k_n0 * exp(p.E.kn/p.R*(1/p.T_ref - 1/TEMP)); 
 p.k_p = p.k_p0 * exp(p.E.kp/p.R*(1/p.T_ref - 1/TEMP)); 

% Exchange current density [I m-2]

 i_0n = p.k_n * p.Faraday * sqrt(((p.c_s_n_max - c_ss_n) .* c_ss_n .* p.ce));
 i_0p = p.k_p * p.Faraday * sqrt(((p.c_s_p_max - c_ss_p) .* c_ss_p .* p.ce));

 xn= 0.5 * (cur/p.Area_n)/(p.a_n*p.L_n)/i_0n;
 xp=-0.5 * (cur/p.Area_p)/(p.a_p*p.L_p)/i_0p;

% Overpotentials [V]
 RTaF=(p.R*TEMP)./(p.alph*p.Faraday);
 
%  eta_n  = RTaF .* asinh(cur ./ (2.*p.a_n.*p.Area_n.*p.L_n.*i_0n));
%  eta_p  = RTaF .* asinh(-cur ./ (2.*p.a_p.*p.Area_p.*p.L_p.*i_0p));

 eta_p=(2*p.R*TEMP)/p.Faraday * log(xp + sqrt(1+xp*xp));
 eta_n=(2*p.R*TEMP)/p.Faraday * log(xn + sqrt(1+xn*xn));

% SPM Voltage [V]
 V_spm= eta_p - eta_n + V_ocv;


%% Degredation  


Sn=3*p.eps_s_n*p.L_n*p.Area_n/p.R_n;  %m^2                                        
                              
% Overpotential of the SEI [V]

eta_sei_n= Uref_n + eta_n - p.Us + (sei*2037.4)*cur;                                 

% Current density of the SEI [I m-2]

ksei = p.ksei* exp((130000/p.R)*(1/p.T_ref - 1/TEMP)); 
is= real (p.Faraday*ksei*exp(-p.Faraday/(p.R*TEMP) *( eta_sei_n) ));                           

% Growth rate of SEI layer [m]

sei_dot = is/(p.Faraday*p.rhos);

% Total resistance (film + growing SEI layer)

R_tot_n = p.Rsei_n*sei + p.Rcell;

% cyclable Li, evolution of the charge [C]

Q_dot= -Sn*is;                                                                   

% Total molar flux [mol m-2 s-1]

js_n = (J_n) + (is/p.Faraday);

% Lithium Loss [A]

LLi=is*p.Area_n*p.L_n*p.a_n;

% Solid particle concentration [mol m-3]

c_p = A_p*U_p + B_p.*J_p;
c_n = A_n*U_n + B_n.*js_n;  


%% Heat generation  


% Heat remove
 [Qohmic,Qremv]=coolingcontrol(T,cur,R_tot_n,p);
 Qgen= Qohmic + 0 + 0;

 % Temperature calculation [K]
 T1_dot= Qgen./p.Cc + (T2-T1)./(p.Rc*p.Cc); %Tc core temp
 T2_dot= (p.T_cool-T2)./(p.Ru*p.Cs) - (T2-T1)./(p.Rc*p.Cs); %Ts surface tem. 

% Lumped Temperature calculation [K]
 T_dot= (Qgen  - Qremv)./(p.M*p.Cp);
 
 
%% Outputs
xdot = [c_n; c_p; real(T1_dot); real(T2_dot); real(T_dot); real(Q_dot); real(sei_dot)]; 


varargout{1} = theta_p;
varargout{2} = theta_n;
varargout{3} = V_spm;
varargout{4} = V_ocv;
varargout{5} = eta_n;
varargout{6} = LLİ;
varargout{7}= Uref_n;
varargout{8}= Qohmic;
varargout{9}= Qremv;
varargout{10}= R_tot_n;
varargout{11}= cur;
end


function [check,isterminal,direction]=Efcn(t,x)
global p

check= x(49) <= p.c_s_n_max*p.theta_n_min || x(98) >= p.c_s_p_max*p.theta_p_max ; 
direction=[];
isterminal=1;
    
end

function [check,isterminal,direction]=Efcn1(t,x)
global p

check= x(49) >= p.c_s_n_max*p.theta_n_max || x(98)<= p.c_s_p_max*p.theta_p_min;
direction=[];
isterminal=1;

end
