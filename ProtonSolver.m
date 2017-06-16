%this is the velocity / LET solver approximation, version 2: uses EiInv
%clear 
close all 

disp('Proton Energy Bethe Equation solver: David Robert Grimes and Daniel Warren, University of Oxford.')
disp('Please note this code is provided as a simple proof of concept, and is not optimised!')
disp('E-mail: davidrobert.grimes@oncology.ox.ac.uk / grimesd2@gmail.com Twitter: @drg1985') 
%Constants
me = 9.10938356e-31; %electron mass
e = 1.60217662e-19; %electron charge
eo = 8.854187817e-12; %vacuum permittivity 
I = 75.*e; % ionization mean
n = 3.343e29; %electron density (electrons/cc)
mp = 1.6726219e-27; %particle mass - in this case proton
z = 1; %charge in multiples of electron
c = 299792458; %speed of light

pe = input('Enter Initial Particle Energy in MeV : ');
ke_p = pe*1e6*e; % proton energy (converted from MeV ie. 150 MeV)
gam_p = (mp*c^2+ke_p)/(mp*c^2); % relativstic gamma
bet_p = sqrt(1-(1/gam_p)^2); % relativistic beta

vo = bet_p.*c ; %initial particle velocity

%Collect together unweldy terms! 

fc = (e^2)./(4*pi*eo);
A = ((4*pi*n*(z^2))./me).*(fc^2);
B = (2*me)/I;
As = A./mp; %ignore minus



load('EiInvTotalFull.mat');



% START OF ANALYTICAL SOLUTION

d = 0:0.0001:1; % d is depth in matter in units metres, millimetres here

R = ei(log(B^2*vo^4))/(2*B^2*As); % calculate max range

h = waitbar(0,'Finding inverse exp. ints...');
for i = 1:numel(d)
    
     disE = R-d(i);
     TablePlace = round(2*B^2*As*disE); %Find position in look up table  
    
    if TablePlace > 0
        
         
        
           
        inv = EiInv(TablePlace);
        v(i) = ((1/B^2)*exp(inv))^.25;
        h = waitbar(i/numel(d),h);
    else
        v(i) = 0; % speed is zero beyond end of range
         h = waitbar(i/numel(d),h);
    end
end
close(h);

%now we make a gamma cubed function


gamma = 1./sqrt(1 - (v./c).^2);
gamma3 = gamma.^3;
g3 = spline(d,gamma3);

%transform x' -> x 

u = length(d); 
for y = 1:u
xt(y) = integral(@(x)ppval(g3,x),0,d(y)); 
end

Rn = integral(@(x)ppval(g3,x),0,R);
MaxRPlot = 1.1*Rn; 

%figure

%Energy of particle

gn = 1./sqrt(1 - (v./c).^2);
Energy = (gn-1).*mp.*(c^2); 
EnergyMEV = Energy./(10^6*e);


%Particle LET 

LET = (A./(v.^2)).*log(B.*v.^2);

LET = (10^-9).*(LET./e); %LET in Kev / um; 

%Create plots of energy, velocity and LET
str1 = 'Energy Profile for ';
str2 = 'Velocity Profile for ';
str3 = 'LET Profile for ';
strid = ' MeV proton';

vfrac = v./c; 



subplot(2,2,[1 2]); plot(xt,EnergyMEV); xlabel('Distance in tissue (m)'); ylabel('Particle Energy (MeV)'); title([str1 num2str(pe) strid]); xlim([0 MaxRPlot]);
subplot(2,2,3); plot(xt,vfrac); xlabel('Distance in tissue (m)'); ylabel('Particle velocity (Fraction of light speed)'); title([str2 num2str(pe) strid]); xlim([0 MaxRPlot]);
subplot(2,2,4); plot(xt,LET); xlabel('Distance in tissue (m)'); ylabel('LET (keV / um)'); title([str3 num2str(pe) strid]);xlim([0 MaxRPlot]);
 set(gcf, 'Position', get(0, 'Screensize'));
%clearvars -except LET xt EnergyMEV v Rn 
