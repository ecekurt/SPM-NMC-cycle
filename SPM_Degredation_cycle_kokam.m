%% KOKAM NMC battery SPM model with SEI effect
% Fitted electrode thicknesses are used. Other parameters are same with the
% SLIDE Kokam NMC parameters

clear all
close all
clc 

run Kokam_NMC_parameters
load('KokamOCVNMC.mat'); %half cell pot.
load('KokamOCVC.mat');  %half cell pt.
load('OCVcell.mat');  %overall cell ocv vs Ah
load('KokamNMC.mat'); 
load('KokamC.mat');   
global data 
global p
global KokamOCVNMC KokamNMC 
global KokamOCVC KokamC
global OCVcell

%% Constant current input

% Discharge only
 p.C_rate= 1;
 p.cycle_number=5;
 p.Nc= p.cycle_number*2;

tend= 3600*p.Nc/p.C_rate;
data.Current=@(t) 0*(t==0) - 2.7*(t<=1) + p.C_rate*2.7*(1<=tend);

%% Finite difference for spherical particle and electrolyte
p.Np=50;
p.Nn=50;
p.delta_p =  p.R_p/(p.Np);
p.delta_n =  p.R_n/(p.Nn);  

%% Initial concentration  of solid particles and electrolyte 


Up0=p.c_s_p_max*p.theta_p_min*ones(p.Np-1,1);             
Un0=p.c_s_n_max*p.theta_n_max*ones(p.Nn-1,1);            
p.T_amb=298.15;
% Temperature
T10 = p.T_amb;
T20 = p.T_amb;
T0 =  p.T_amb;

% SEI

Qs0=p.eps_s_n*p.Faraday*p.Area_n*p.L_n*p.c_s_n_max*p.theta_n_max;
sei0=p.L_sei;

 tic   
% set up starting equation    
tspan=0:1:tend;
t=tspan(1);
x = [Un0; Up0; T10; T20; T0;Qs0; sei0]';
func=@ode_SPMT_discharge;
options=odeset('Events',@Efcn); 

eventtime=[];
eventtime(1)=0;
a=0;
    
while t(end) < tend % main loop
      
if a>=p.cycle_number*2
    break
end


% Run integration until event function stops it
[at,ax,ate,aye,aie] = ode23s(func,[length(t):tend],x(end, :),options); 

if ~isempty(ate)
eventtime=[eventtime,round(ate)];
end

% Append the new trajectory
t = cat(1, t, at(2:end)); 
x = cat(1, x, ax(2:end,:)); 


% Decide for new function and event function
if (x(end,49) <= p.c_s_n_max*p.theta_n_min || x(end,98) >= p.c_s_p_max*p.theta_p_max )
    func=@ode_SPMT_charge;
    options=odeset('Events',@Efcn1);    % N>=25000 || P<=680  buradan 680den küçük olunca çıkacak dolayısıyla cn her türlü 1000den büyük olacak
    
    elseif(x(end,49) > 1000)
    func=@ode_SPMT_discharge;
    options=odeset('Events',@Efcn);     % N<=25 || P>=18000  
    
end

 a=a+1;
   
end




c_n = x(:,1:(p.Nn-1));
c_p = x(:,p.Nn : 2*(p.Nn-1));
T1 = x(:,end-4);
T2 = x(:,end-3);
T  = x(:,end-2);
Capacityloss= x(:,end-1);
seigrowth= x(:,end);

% Extract data
for i=1:p.Nc
  duration=diff(eventtime);
if mod(i,2)==1
        
    
    for k=1+eventtime(i):eventtime(i+1)

    [~,theta_p(k),theta_n(k),V_spm(k),V_ocv(k), Uref_n(k), Uref_p(k), ...
        cur(k)]...
        =ode_SPMT_discharge(t(k),x(k,:)');
    SOCp(k)=( theta_p(k)- p.theta_p_max )/( p.theta_p_min -p.theta_p_max);
    SOCn(k)=( theta_n(k)- p.theta_n_min )/( p.theta_n_max -p.theta_n_min);
    
    end

  
else
      
        
    for k=1+eventtime(i):eventtime(i+1)
        
    [~,theta_p(k),theta_n(k),V_spm(k),V_ocv(k), Uref_n(k), Uref_p(k), ...
        cur(k)]...
        =ode_SPMT_charge(t(k),x(k,:)');
    SOCp(k)=( theta_p(k)- p.theta_p_max )/( p.theta_p_min -p.theta_p_max);
    SOCn(k)=( theta_n(k)- p.theta_n_min )/( p.theta_n_max -p.theta_n_min);
    
    end
       
end
Ah(i)=integral(data.Current,0,duration(i));

end

Qp=(p.Area_p*p.L_p*p.Faraday*p.eps_s_p*p.c_s_p_max*p.theta_p_max)/3600;
Qn=(p.Area_n*p.L_n*p.Faraday*p.eps_s_n*p.c_s_n_max*p.theta_n_max)/3600;


Ah_th=cumsum(Ah)/3600;
Nc=linspace(0,Ah_th(end),length(Capacityloss));
Qloss=(Capacityloss/3600);
plot(Nc,(100*Qloss)/Qn,'.-'); xlabel('Ah - processed [Ah]');ylabel('Charge loss [%]');
grid on;

cycle=linspace(0,p.cycle_number,length(V_spm));
toc

figure
plot(cycle,V_spm,'.-');
xlabel('Cycle Number');
ylabel('Voltage [V]');
legend('');
grid on;

figure
plot(cycle,theta_n,'.-');
xlabel('Cycle Number');
ylabel('Negative electrode Li+ fraction');
grid on;

figure
plot(cycle,seigrowth,'.-'); 
xlabel('Cycle Number');
ylabel('Sei thickness [m]');
grid on;


function [xdot,varargout]=ode_SPMT_discharge(t,x)

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


cur=data.Current(t);


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


%% Degredation  


Sn=3*p.eps_s_n*p.L_n*p.Area_n/p.R_n;  %m^2                                      

% Overpotential of the SEI [V]

eta_sei_n= Uref_n + eta_n - p.Us + (sei*p.Rsei_n)*cur;       %2037.4 R                        

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

% Solid particle concentration [mol m-3] 

c_p = A_p*U_p + B_p.*J_p;
c_n = A_n*U_n + B_n.*js_n;  

%% Heat generation  

Qohmic =  cur.^2*R_tot_n;
Qgen= Qohmic + 0 + 0;
  
% Heat remove
Qremv= p.h*p.A*(T-p.T_cool);  

% Temperature calculation [K]
T1_dot= Qgen./p.Cc + (T2-T1)./(p.Rc*p.Cc); %Tc core temp
T2_dot= (p.T_amb-T2)./(p.Ru*p.Cs) - (T2-T1)./(p.Rc*p.Cs); %Ts surface tem. 

% Lumped Temperature calculation [K]
T_dot= (Qgen  - Qremv)./(p.M*p.Cp);

%% Outputs
xdot = [c_n; c_p; real(T1_dot); real(T2_dot); real(T_dot); real(Q_dot); real(sei_dot)]; 

varargout{1} = theta_p;
varargout{2} = theta_n;
varargout{3} = V_spm;
varargout{4} = V_ocv;
varargout{5} = Uref_n;
varargout{6} = Uref_p;
varargout{7}= cur;
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


cur=-data.Current(t);


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

% li-fraction  of the electrodes 
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
 
 eta_n  = RTaF .* asinh(cur ./ (2.*p.a_n.*p.Area_n.*p.L_n.*i_0n));
 eta_p  = RTaF .* asinh(-cur ./ (2.*p.a_p.*p.Area_p.*p.L_p.*i_0p));

%  eta_p=(2*p.R*TEMP)/p.Faraday * log(xp + sqrt(1+xp*xp));
%  eta_n=(2*p.R*TEMP)/p.Faraday * log(xn + sqrt(1+xn*xn));

% SPM Voltage [V]
 V_spm= eta_p - eta_n + V_ocv;


%% Degredation  


Sn=3*p.eps_s_n*p.L_n*p.Area_n/p.R_n;  %m^2                                                                         

% Overpotential of the SEI [V]

eta_sei_n= Uref_n + eta_n - p.Us + (sei*p.Rsei_n )*cur;                                 

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

% Solid particle concentration [mol m-3]

c_p = A_p*U_p + B_p.*J_p;
c_n = A_n*U_n + B_n.*js_n;  


%% Heat generation  

Qohmic =  cur.^2*R_tot_n;
Qgen= Qohmic + 0 + 0;
  
% Heat remove
Qremv= p.h*p.A*(T-p.T_cool);  

% Temperature calculation [K]
T1_dot= Qgen./p.Cc + (T2-T1)./(p.Rc*p.Cc); %Tc core temp
T2_dot= (p.T_amb-T2)./(p.Ru*p.Cs) - (T2-T1)./(p.Rc*p.Cs); %Ts surface tem. 

% Lumped Temperature calculation [K]
T_dot= (Qgen  - Qremv)./(p.M*p.Cp);

%% Outputs
xdot = [c_n; c_p; real(T1_dot); real(T2_dot); real(T_dot); real(Q_dot); real(sei_dot)]; 

varargout{1} = theta_p;
varargout{2} = theta_n;
varargout{3} = V_spm;
varargout{4} = V_ocv;
varargout{5} = Uref_n;
varargout{6} = Uref_p;
varargout{7}= cur;
end





function [check,isterminal,direction]=Efcn(t,x)
global p

check= x(49) <= p.c_s_n_max*p.theta_n_min || x(98) >= p.c_s_p_max*p.theta_p_max ; %theta_n_min=0.0316 Vmin=2.7 
direction=[];
isterminal=1;
    
end

function [check,isterminal,direction]=Efcn1(t,x)
global p

check= x(49) >= p.c_s_n_max*p.theta_n_max || x(98)<= p.c_s_p_max*p.theta_p_min;
direction=[];
isterminal=1;

end
