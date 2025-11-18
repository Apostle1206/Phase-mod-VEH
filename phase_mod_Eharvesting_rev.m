function phase_mod_Eharvesting

%Data files
%Eh_QvsSOM_FORCE_Symm_sawtooth_g1pt5.dat
%Eh_QvsSOM_FORCE_square_wave_g1pt2.dat
%not used Eh_QvsSOM_FORCE_asymmetric_sawtooth_g1pt5
%Eh_QvsSOM_FORCE_cor_asymmetric_sawtooth_g1pt5.dat
%Eh_QvsSOM_FORCE_periodic_cosine_wave_g1pt5
fid=fopen('Eh_QvsSOM_FORCE_square_wave_g1pt2.dat','w');
global h;
tic;
Kc= 4.15e3;
Kv=0.0011;
g=1.2;
V0=1;
lambda=0.5;
f=0.1;
R=300e3;
Cp=112e-9;
Tp=R*Cp;
d=2.97;
mu=pi*4e-7;
mm= 0.051;
A=d^2*((mu*mm^2)./(2*pi*d))^(-2/3);
B=A./(d^2);

somList=0.01:0.01:4;

n1=10.;n2=100.;
Y0 =[0.0 0.0 0.0];
opts = odeset('RelTol',1e-8,'AbsTol',1e-10 );
for j=1:length(somList)
    som=somList(j); 
    
    bom=6.7*som;
    period=2.*pi/som;
step=2.*pi/som/100.;

[~, Y] = ode45(@(t, y) vrf(t, y, g, bom), [0 n1 * period], Y0, opts);
    Y1 = Y(end, :);
   [t, Y] = ode45(@(t, y) vrf(t, y, g, bom), 0:step:n2 * period, Y1, opts);
    
        

  vol=Y(:,3).*Y(:,3)/R;      
z1=trapz(t,Y(:,1).*cos(som*t),1);
z1=2.*z1/n2/period;
z2=trapz(t,Y(:,1).*sin(som*t),1);
z2=2.*z2/n2/period;
q1=(sqrt(z1*z1+z2*z2))/f;

s1=trapz(t,vol.*cos(som*t),1);
s1=2.*s1/n2/period;
s2=trapz(t,vol.*sin(som*t),1);
s2=2.*s2/n2/period;
q2=(sqrt(s1*s1+s2*s2))/f;
AvgP=mean(vol);
%U=[U;q];
%Power output Px
xt1=trapz(t,Y(:,1),1);
xt1=2.*xt1/n2/period;
Ux=xt1;
P_outx=trapz(t,(Y(:,1)-Ux).^2 ,1);
P_outx=2.*P_outx/n2/period;

%Power output Pv
P_outv=trapz(t,(vol).^2 ,1);
P_outv=2.*P_outv/n2/period;

%Power input
P_in=(f^2)/2;

%Efficiency in x
effx=(P_outx)/(P_in);

%Efficiency in v
effv=(P_outv)/(P_in);


fprintf(fid,'%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f  %12.8f\n', som,q1,q2, AvgP,P_outx,P_outv, P_in,effx, effv);
end
%fprintf(fid,'\n');
toc;


%-----------
function dY = vrf(t,V,g,bom)
 %  The derivative function
x=V(1);y=V(2); z=V(3);
dx=y; 


% Call the asymmetric_sawtooth function
%b = newest_symmetric_sawtooth(t, bom, g);

% Call the square_wave function
h = new_square_wave(t, bom, g);

% Call the asymmetric_sawtooth function
%P = cor_asymmetric_sawtooth(t, bom, g);

%Define the periodic cosine function
%rr = periodic_cosine_wave(t, bom, g);

dz=Kc*y - z/Tp;

dy=-lambda*y-V0*sin(x)*cos(h)-V0*cos(x)*sin(h)...
    - Kv*z +f*cos(som*t);
dY=[dx;dy;dz];

end
end