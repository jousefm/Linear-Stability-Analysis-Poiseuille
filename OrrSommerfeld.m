clear all; 

close all; 

clc

%% =======================================
% Jousef Murad 3rd Semester - KIT
% Numerical Fluid Mechanics II
% Supervisor: Prof. Dr. Markus Uhlmann
% ========================================
% Solution of the Orr-Sommerfeld equation
% ========================================

%% Flow Parameter
tic

a = 1.0; %+h

U0 = 1.0; % Ref. Velocity / Max. Velocity at y = 0

%% Discretzation

N = 200; 

%% Define Discretization

Dy = 2/N; % Grid discretization taking into account "2h" which is the height of the channel

imUnit = sqrt(-1); % not i because of index in loops


%% Velocity and Derivative

% Pre-Allocation for y probably needed --> test how much faster it is
% without Pre-Allocation

for i = 1 : N+1
   
    y(i) = (i-1.0)*Dy - 1.0;
    
    U(i) = U0*(1.0-(y(i))^2);
    
    UDiff(i) = -2.0*U0;   % 2nd Derivative
    
end

%% Matrices for 4th and 2nd derivatives from OS

%% Matrix for 4th derivative f''''
% Pentadiagonal (Abbreviation: PD)

PD = zeros(N-1,N-1); % Pre-Allocation of matrix

PD(1,1) = 7; % will be discussed in the small thesis why this is 7 --> Related to BCs
PD(1,2) = -4; 
PD(1,3) = 1;
PD(2,1) = -4; 
PD(2,2) = 6; 
PD(2,3) = -4; 
PD(2,4) = 1;

for i = 3 : N-3
    
    PD(i,i-2) = 1; PD(i,i-1) = -4; PD(i,i) = 6; PD(i,i+1) = -4; PD(i,i+2) = 1;
    
end

PD(N-2,N-4)=  1.0; 
PD(N-2,N-3)= -4.0; 
PD(N-2,N-2)=  6.0;
PD(N-2,N-1)= -4.0;
PD(N-1,N-3)=  1.0; 
PD(N-1,N-2)= -4.0; 
PD(N-1,N-1)=  7.0;
        
PD = PD/Dy^4; %IMPORTANT!!!


%% Matrix for 2nd derivative f''
% Tridiagonal (Abbreviation: TD)

TD = zeros(N-1,N-1); % Pre-Allocation of matrix

TD(1,1) = -2.0; 
TD(1,2) = 1.0;

for i=2:N-2
    
  TD(i,i-1)=1.0; TD(i,i)=-2.0; TD(i,i+1)=1.0;
  
end

TD(N-1,N-2) = 1.0; TD(N-1,N-1) = -2.0;

TD = TD/Dy^2; %IMPORTANT!!!

%% Matrix for 2nd derivative U*f''
% Tridiagonal (Abbreviation: T1)

% Error from my side --> Initially set U(1) as start value and even forgot
% to multiply with the velocity so that the plots looked horrible and
% started at 1 instead of 0.65

T1 = zeros(N-1,N-1); % Pre-Allocation of matrix

T1(1,1) = -2*U(2); %U(1) = 0
T1(1,2) =  1*U(2);

for i=2:N-2
  
  UStar = U(i+1); %U(2) already covered top and need to incrementally step forward in the velocity profile
    
  T1(i,i-1) = 1.0*UStar; T1(i,i) = -2.0*UStar; T1(i,i+1) = 1.0*UStar;
  
end

T1(N-1,N-2) =  1.0*U(N); 
T1(N-1,N-1) = -2.0*U(N);

T1 = T1/Dy^2; %IMPORTANT!!!

%% Wave number definition and evalutation

kIter = 50; % Wave numbers for the loop

%k = 1.0
%Deltak = k/kIter;

looptoken = 0;

if (looptoken == 0)

ReArray = [6000 , 3000 , 6000];
alpha = [1 , 0.5 , 1 , 2, 5];

else
  
alpha = 0:1:32;
ReArray = 0:1:32;  
    
end

numRe = 50; % Amount of Re numbers scanned during the process

ReStart = 2000;

ReEnd = 40000;

Dalpha = 1.2/10; % Maximum wave number investigated over iteration of loops of alpha/k

numReyloop = 10; % Also determines the "fineness" of the 2D contourplot like denominator of Dalpha

%Ramp = 1.5; % Ramping for the Reynolds Number

Ramp = (ReEnd/ReStart)^(1/(numReyloop-1));

%% FOR-LOOP

for idx = 1:11 % loop over wave number

%alpha1 = alpha(idx);

alpha1 = Dalpha*(idx - 1.0) + 0.0001;

%% For loop for Reynolds number

%Re = ReArray(1);

Re = 2000; % Set for kin. viscosity

%ReIdx = 1;

for idx2 = 1:10 % loop over Reynolds number
    
%ReNew = ReArray(idx2); 

% if (idx2 == 1)
%     
% nue = 1/Re; % kinematic viscosity in [m^2/s]
% 
% else 
%     
%     nue = 1/Re(idx2)
%     
% end

nue = 1/Re;

A = - 2.0*alpha1*imUnit*nue*TD + imUnit*nue/alpha1 * PD + T1;
   
   
for i=1:N-1
    
    A(i,i) = A(i,i) - UDiff(i+1) - U(i+1)*alpha1^2 + imUnit*nue*alpha1^3;
    
end

B = TD;

for i=1:N-1
    
    B(i,i) = B(i,i)-(alpha1*alpha1);
    
end

%% Eigenvalues from system

[V,D] = eig(A,B); 
% Returns diagonal matrix D of generalized eigenvalues and full matrix V whose 
% columns are the corresponding right eigenvectors, so that A*V = B*V*D

egvl = eig(A,B);
% Returns a colum vector containing the generalized eigenvalues of square
% matrices A and B

G2 = imag(alpha1*V); % Growth rate

G = imag(alpha1*egvl);

GSort = -sort(-G); % Sort Growth Rate - Last value biggest value !!!

GMAX(idx,idx2) = GSort(1);

wavenum(idx) = alpha1*a; %Wavenumber

Re = Ramp * Re;

%% Eigenvalue Spectrum

egvlspec = 0;

if (egvlspec ==1)
    
for j=1:N-1;
    
realc(j)=real(D(j,j));
imagc(j)=imag(D(j,j));
    
end

figure(1)
hold on
plot(realc,imagc,'+',[0 1],[0 0])
axis([0 1 -1 0.1]) % set x-max to 1 --> still wrong...
%legend('k = 0.1' , 'k= 0.5' , 'k = 1' , 'k = 2', 'k = 5')
title('Eigenvalue Spectrum - Re = 6,000 - Alpha = 1','FontSize',25)
xlabel('c_r', 'FontSize', 25)
ylabel('c_i', 'FontSize', 25)
set(gca,'fontsize',25)

end

Rey(idx2) = Re/Ramp; %Used for iteration of the stability contour plot

end



end


% ====================
% Plot Section =======
% ====================

%% Plot Tokens

loglogplot = 0;
plot2dtoken  = 0;
contplottoken = 1;
timetoken = 0;
accuracytoken = 0;

i = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(loglogplot ==1)
    
Discret = [100 200 300 400 500];
L2errImag = [0.003125654888904 0.003640397573013 0.003699963536898 0.003718194299598 0.003726178966810];
Literature = 0.00373967; %Literature Value of the Imaginary Part for Alpha = 1, Re = 10,000 according to Orszag using Chebyshev polynomials 
                         %Compare with GSort variable
                         
relError = (1 + (Literature-L2errImag./Literature))*100;

figure(i+1)
hold on

p = polyfit(log(Discret), log(relError), 1);

set(gca, 'XScale', 'log','YScale', 'log')
loglog(Discret,relError,'LineWidth',3)
            
triang_x = [Discret(2), Discret(3)];  %Location of triangle start
triang_y = interp1(Discret, relError, triang_x);
loglog(triang_x([1,2,2]), triang_y([1,1,2]), 'k', 'LineWidth', 2);

x1 = 320; %Vector for annotation
y1 = 2; %Vector for annotation

str1 = sprintf('m = -1.96');

text(x1(1),y1(1),str1,'FontSize', 25)

grid on

title('Convergence Study','FontSize',25)
xlabel('Discretization', 'FontSize', 25)
ylabel('Relative Error', 'FontSize', 25)
set(gca,'fontsize',25)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(plot2dtoken==1)
    
figure(i+2)
hold on
set(gca,'fontsize',15)
box on
xlabel('ka','fontsize',15)
ylabel('\sigma_I{a}/U_0','fontsize',15)
% axis([0 1.2 -0.06 0.01])

plot(wavenum,GMAX(:,1),'k','linewidth',3);

for idx2=2:ReIdx
    plot(wavenum,GMAX(:,idx2),'k');
end

plot([0 1.2],[0,0],'--k')

end

%% Plot 2

if(contplottoken==1)
    
figure(i+3)
hold on
set(gca,'fontsize',25)
box on
xlabel('Re','fontsize',25)
[C,h]=contour(Rey,wavenum,GMAX,[0.0, 0.0],'k-');
plot([5850 5850],[0.55,1.2],'k--')
ylabel('ka','fontsize',25)

x1 = [25000 18000]; %Vector for annotation
y1 = [1.1 0.88]; %Vector for annotation

str1 = sprintf('stable');
str2 = sprintf('unstable')

text(x1(1),y1(1),str1,'FontSize', 25)
text(x1(2),y1(2),str2,'FontSize', 25)

axis([0 40000 0.55 1.2])
end
%% Time Plot

if(timetoken ==1)
    
Points = [50 100 200 300 400 500];
Time = [0.161663 0.227687 0.3045 0.541223 1.031586 1.887334];

figure(i+5)
plot(Points,Time,'-.k*','LineWidth',3,'MarkerSize',5)
grid on
title('Time to solve Eigenvalue Problem - Re = 10,000 - Alpha = 1','FontSize',25)
xlabel('Discretization','fontsize',25)
ylabel('Time','fontsize',25)
set(gca,'fontsize',25)
axis([50 500 0 2])
end

if (accuracytoken == 1)
    
Discret = [100 200 300 400 500];
L2errImag = [0.003125654888904 0.003640397573013 0.003699963536898 0.003718194299598 0.003726178966810];
Literature = 0.00373967; %Literature Value of the Imaginary Part for Alpha = 1, Re = 10,000 according to Orszag using Chebyshev polynomials 
                         %Compare with GSort variable

relError = (1 + (Literature-L2errImag./Literature))*100;
   
figure(i+5)
plot(Discret,relError,'-.r*','LineWidth',3,'MarkerSize',5)
grid on
title('Accuracy of the imaginary part of most unstable EV - Re = 10,000 - Alpha = 1','FontSize',25)
xlabel('Discretization','fontsize',25)
ylabel('Error [%]','fontsize',25)
set(gca,'fontsize',25)
%axis([0 500 0 2])
end

toc
