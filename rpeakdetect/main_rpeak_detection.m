clc;
close all;

record='108';

sfile=strcat('msmecg\msm',num2str(record),'.dat');

x=load(sfile);
% x=awgn(x,30);
x1=x(:,1);

% x1= medfilt1(x1,10);

bm1=0;

ND=floor(20*360);

MF=ceil(length(x1)/ND);

N1=bm1*ND+1;

for bm=bm1+1:MF-1

N2=bm*ND;


Fs=360; 
N=N2-N1;

x=x1(N1:N2);

%-----------------------------------------

FR=360;

nL=length(x);

[see, d, rlap]=mfun_peakdetection(x,Fs);

%-----------------------------------------

xs=x./max(abs(x));
x_m1=xs;
s2=smooth(x_m1,200,'moving');
sig1=x_m1-s2;

sig1=sig1./max(abs(sig1));

QRS_S=rlap;

NW=30;

NW1=30; 

% ws=2*NW1+1;
% xn(1:2:ws)=1;
% xn(2:2:ws)=-1;

nq=1; 

sig2=[zeros(1,NW1)'; sig1; zeros(1,NW1)'];

QRS_S=QRS_S+NW;

for k=1:length(QRS_S)
            
     WI=QRS_S(k)-NW;
     
     WE=QRS_S(k)+NW;
     
     wsig=sig2(WI:WE);
     
     [R_VS1, R_L1]=max(abs(wsig));
     
%      [R_V1 RL2]=max(see(WI:WE));
     
%      wsig2=sig1(R_L1+WI-NW1:R_L1+WI+NW1);
     
       
     R_W(nq,:)=[R_VS1  R_L1+WI-1 ];
     
     nq=nq+1;
     
end

B_L=R_W(:,2)-NW;

annotf=strcat('QRS_L\',num2str(record),'.mat');
load(annotf)
 
A_L=L;
f=find(N1 < A_L & A_L < N2);
A_L1=A_L(f);
A_L1=A_L1-N1+1;

t=(0:length(sig1)-1)./FR;
rp1=ones(1,length(rlap));

see=see(1:ND);

d=d(1:ND);


bm;


 figure;subplot(311);plot(x);axis tight;grid on;
 subplot(312);plot(d);axis tight;grid on; ylabel('r3')
 subplot(313);plot(see);axis tight;grid on;ylabel('nzcrp')



figure(3);subplot(511);plot(t,x);axis tight;grid on;ylabel('x[n]');xlabel('(a)');
subplot(512);plot(t,d);axis tight;grid on;ylabel('d[n]');xlabel('(b)');
subplot(513);plot(t,see);axis tight;grid on; ylabel('s[n]');xlabel('(c)');
subplot(514);stem(t(B_L), x(B_L),'r'); hold on; plot(t,x);axis tight; grid on; ylabel('Detected R-peaks');
subplot(515);stem(t(A_L1), x(A_L1),'r'); hold on; plot(t,x);axis tight; grid on; ylabel('Detected R-peaks');

xlabel('Time (sec)');

f1 = figure(3);
scrsz = get(0,'ScreenSize');
set(f1, 'Position', [1  1  scrsz(3) scrsz(4) ] );




pause;close all; 

clear B_L A_L B_L1 QRS_S z1 R_W 

N1=N2+1;

end

