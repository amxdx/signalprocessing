function red1_index=peakdetection(p)
%sample rate of the signals=360Hz
Fs=360;
n=load(p); %1min record containing normal beat waveforms
%above ecg .mat files contains 2 channel ECG records.
%we only need one channel so extracting only one of the channels
%(channel MLII)
normal= n.val(1,:);
%normal=normal/max(abs(normal));

%%%% plots %%%%

plot(normal);
title('Raw ECG Training Data (normal)')
normal0=normal;
%% Filtering%%
% % Removing DC Components of the ECG Signal
% normal=normal-mean(normal);
% 
% % Removing  High Frequency Noise
% B1=(1/10)*ones(1,10);
% A1=1;
% normal=filter(B1,A1,normal);
% 
% %Removing  High Frequency Noise
% B2=(1/1.0025)*[1 -1];
% A2=[1 -0.995];
% normal=filter(B2,A2,normal);
% 
% % Removing 60Hz Powerline Interference(Comb Filter)
% B3=conv([1 1],[0.6310 -0.2149 0.1512 -0.1288 0.1227 -0.1288 0.1512 -0.2149 0.6310]);
% A3=1;
% normal=filter(B3,A3,normal);

% TV
Nit=20;
lambda=2;
y=normal/max(abs(normal));
J = zeros(1,Nit);
[ytv,J] = denoiseTV(y,lambda,Nit);
figure;plot(ytv);
title('ECG Data after preprocessing');
% normal1=normal;

 
%% 
ytv_diff=diff(ytv);
ytv_enr=ytv_diff.^2;
% thresh0=2*(std(ytv_enr));
% ytv_enr=(ytv_enr>thresh0).*(ytv_enr);
Nb=floor(0.1*Fs);
b=rectwin(Nb)/Nb;

a=1;
ytv_enr=filtfilt(b,a,ytv_enr);
x=ytv_enr/max(ytv_enr);
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
figure;
plot(normal);
for i = 1:h
    hold on;
    if x1(i)~=0
  plot(i,x1(i),'o','color','red')
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
subplot(2,1,1);plot(x);    
subplot(2,1,2);plot(x); hold on 
subplot(2,1,2);plot(index1,x1(index1),'+','color','green');
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
        


figure;    
subplot(2,1,1);plot(x);    
subplot(2,1,2);plot(x); hold on 
subplot(2,1,2);plot(red1_index,x1(red1_index),'+','color','red');
 
thresh1=1/5*(mean(x1(red1_index)));
for i=1:length(red1_index)
if x1(red1_index(i)) < thresh1
    x1(red1_index(i))=0;
end
end
for i=1:length(red1_index);
indx_g1=find(x1(red1_index));
end

figure;    
subplot(2,1,1);plot(x);    
subplot(2,1,2);plot(x); hold on 
subplot(2,1,2);plot(red1_index(indx_g1),x1(red1_index(indx_g1)),'+','color','black');

end