function [x d rpeak]=mfun_peakdetection(sig,Fs)
normal=sig;
% TV
Nit=20;
lambda=2;
y=normal/max(abs(normal));
y=y';
J = zeros(1,Nit);
[ytv,J] = denoiseTV(y,lambda,Nit);
figure;plot(ytv);
title('ECG Data after preprocessing');
% normal1=normal;
% %% 
% d=[0 diff(ytv)];
% ds2=abs(d).^2;
% % ds2=(abs(ds2)>0.06).*ds2;
% ds2=-ds2.*log(ds2);
% WS=45;
% h=rectwin(WS)/WS;  % method I impulse response
% x=filtfilt(h,1,ds2);


d=[0 diff(ytv)];
ytv_enr=d.^2+eps;
% ytv_enr=-ytv_enr.*log(ytv_enr);
% thresh0=2*(std(ytv_enr));
% ytv_enr=(ytv_enr>thresh0).*(ytv_enr);
Nb=floor(0.1*Fs);
b=rectwin(Nb)/Nb;

a=1;
ytv_enr=filtfilt(b,a,ytv_enr);
x=ytv_enr/max(ytv_enr);

%figure;plot(211);plot(x);

%%
max_val = max(x);
tresh = 0.01*max_val;
h = length(x);
p=1;
q=2;
a=0;
for i = 2:h
    if abs(x(p)-x(i))> tresh
        x1(i)=x(i);
        p=p+a+1;
        a=0;
    else
        a=a+1;
    end
end
alp=[];
index2=[];
sgn=[];
sgn(1)=0;
alp(1)=0;
beta=2;
gamma=1;
h=length(x1);

for i = 1:h
    if x1(i)~=0
    alp(beta)=x1(i);
    gamma(beta)=i;
    sgn(beta)=sign(alp(beta)-alp(beta-1));
     if beta > 2
         s_pt= sgn(beta)- sgn(beta-1);
         if abs(s_pt)==2
             x3(gamma(beta))=x1(i);
         end
     end
      beta=beta+1;
    end
end
for i=1:length(x1);
indx_r=find(x1);
end
for i=1:length(x3);
indx_g=find(x3);
end
for i = 1:length(indx_g)
    index(i) = find(indx_r==indx_g(i));
end
index=index-1;
index1=indx_r(index);

% subplot(2,1,1);plot(x);    
% subplot(2,1,2);plot(x); hold on 
% subplot(2,1,2);plot(index1,x1(index1),'+','color','green');

red_index=1;pos1=0;pos2=0;in=1;

for i=1:length(index1)
        temp = x1(index1(i));
        peak_pt = index1(i);
        for k = 1:50
            if peak_pt>55
            temp2(k)= sign(x1(peak_pt)-x1(peak_pt-k));           
            if temp2(k) == -1
                temp2(k)=0;
            end
            else
                temp2(k)=0;
            end
            if peak_pt<(numel(x1)-55)
      temp3(k)= sign(x1(peak_pt)-x1(peak_pt+k));      
           if temp3(k) == -1
                temp3(k)=0;
           end
            end
            
     pos1=nnz(temp2); pos2=nnz(temp3);   
       if pos1>45 && pos2>45
        red_index(i)=index1(i);
       else
            red_index(i)=0;
            end
        
        end
       pos1=0;pos2=0;
end
        ind=find(red_index);
        red1_index=red_index(ind);
        
figure; subplot(2,1,1);plot(sig);    
subplot(2,1,2);plot(x); hold on 
subplot(2,1,2);plot(red1_index,x1(red1_index),'+','color','red');
 

thresh1=1/3.5*(mean(x1(red1_index)));

for i=1:length(red1_index)
if x1(red1_index(i)) < thresh1
    x1(red1_index(i))=0;
end
end

for i=1:length(red1_index);
indx_g1=find(x1(red1_index));
end

% figure;    
% subplot(2,1,1);plot(x); title('final')   
% subplot(2,1,2);plot(x); hold on 
% subplot(2,1,2);plot(red1_index(indx_g1),x1(red1_index(indx_g1)),'b+');

rpeak=red1_index(indx_g1)

% pause
