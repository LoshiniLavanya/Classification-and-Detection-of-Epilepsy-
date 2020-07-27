%%                              SEIZURE DETECTION
%%


clear all;
close all;
clc;

%%                              EEG SIGNAL
%%


s=load('normal2.txt');          %nor = =512      ab   fs  = 173.61
figure;plot(s);
title('EEG SIGNAL SAMPLES')


%%                  TIMA DOMAIN
%%


fs =512;
N=length(s);
ts = 1/fs;
total_rec_time = N / fs;

time_scale = 0:ts:(total_rec_time - ts);
figure;plot(time_scale,s);
title('EEG SIGNAL WITH TIME DOMAIN')

%%                          Sampling frequency
%%

f_hz = freq(N,fs);
display(length(f_hz));

ECG_fft = abs(fft(s));                 % take FFT for convert time domin in to ferq domin
ECG_fft(1) = 0;

figure,
plot(f_hz,ECG_fft);
xlabel('fre (HZ)');
ylabel('ECG fft')
title('ECG FFT ');

%%                              FILTER
%%



[B,A] = butter(3,(1/(fs*2)),'high');                           % use filters to remove noise
filtered_ECG = filter(B,A,s);
[B,A] = butter(3,(30/(fs*2)),'low');
filtered_ECG1 = filter(B,A,filtered_ECG);
figure

subplot(2,1,1);plot(time_scale,s);
hold on

subplot(2,1,2);plot(time_scale,filtered_ECG1,'r');
hold off

title('TOTAL FILTR SIGNAL');



 %%                             DWT
 %%
 
 
waveletFunction = 'db8';
                [C,L] = wavedec(filtered_ECG1,8,waveletFunction);
       
                cD1 = detcoef(C,L,1);
                cD2 = detcoef(C,L,2);
                cD3 = detcoef(C,L,3);
                cD4 = detcoef(C,L,4);
                cD5 = detcoef(C,L,5); %GAMA
                cD6 = detcoef(C,L,6); %BETA
                cD7 = detcoef(C,L,7); %ALPHA
                cD8 = detcoef(C,L,8); %THETA
                cA8 = appcoef(C,L,waveletFunction,8); %DELTA
                D1 = wrcoef('d',C,L,waveletFunction,1);
                D2 = wrcoef('d',C,L,waveletFunction,2);
                D3 = wrcoef('d',C,L,waveletFunction,3);
                D4 = wrcoef('d',C,L,waveletFunction,4);
                                
                [D5,D6,D7,D8] = sub_bands(L,C,waveletFunction);
               figure
                Beta = D5;
                subplot(4,1,1); plot(time_scale, Beta); title('BETA');
                
                
                Alpha = D6;
                subplot(4,1,2); plot(time_scale,Alpha); title('ALPHA'); 
                
                Theta = D7;
                subplot(4,1,3); plot(time_scale,Theta);title('THETA');%%1:1:length(Theta)
                
                 Delta = D8;
                subplot(4,1,4);plot(time_scale,Delta);title('DELTA');
                
                
%%                              SUB BAND SEPERATION
%%
D6 = detrend(D6,0);
xdft3 = fft(D6);
freq3 = 0:N/length(D6):N/2;
xdft3 = xdft3(1:length(D6)/2+1);
ai = round(max(Alpha));
figure;
subplot(411);plot(freq3,abs(xdft3));title('ALPHA');
[~,AI] = max(abs(xdft3));
AI = AI-1;
fprintf('Alpha:Maximum occurs at %d Hz.\n',AI);

                     
D5 = detrend(D5,0);
xdft2 = fft(D5);
freq2 = 0:N/length(D5):N/2;
xdft2 = xdft2(1:length(D5)/2+1);
bi = round(max(Beta));

subplot(412);plot(freq2,abs(xdft2));title('BETA');
[~,BI] = max(abs(xdft2));
BI = BI-1;
fprintf('Beta:Maximum occurs at %d Hz.\n',BI);

 
D7 = detrend(D7,0);
xdft4 = fft(D7);
freq4 = 0:N/length(D7):N/2;
xdft4 = xdft4(1:length(D7)/2+1);
ti  = round(max(Theta));
subplot(413);plot(freq4,abs(xdft4));title('THETA');
[~,TI] = max(abs(xdft4));
TI = TI-1;
fprintf('Theta:Maximum occurs at %d Hz.\n',TI);


D8 = detrend(D8,0);
xdft5 = fft(D8);
freq5 = 0:N/length(D8):N/2;
xdft5 = xdft5(1:length(D8)/2+1);
di = round(max(Delta));
subplot(414);plot(freq3,abs(xdft5));title('DELTA');
[~,DI] = max(abs(xdft5));
DI = DI-1;
fprintf('Delta:Maximum occurs at %d Hz.\n',DI);

if(DI<=5)
    helpdlg('normal')
    fprintf('NOR')


                                                                                                   
 %%                         SVM CLASSIFICATION
 %%
 else
 
DI  =int8(DI);
AI  =int8(AI);
TI = int8(TI);
BI = int8(BI);

QQ = xlsread('Book2.xlsx');
xdata =QQ(1:2,1);

preict = 1;
ictal = 2;
type = [preict;ictal;];
if (13<=TI)
    fprintf('PRE')
%if ((DI>=3)||((6<AI) &&(AI<20)))
%if ((DI<=3)||((6<AI) &&(AI<20))||((12<BI) && (BI<33))|| ((3<TI) && (TI<8)))
result = 1;
elseif ((DI>=6)&&(40<=BI))
    result = 0;
    fprintf('ICT')
end 


[ani_security]=  SVM_train(result,type,xdata);
  
if(ani_security ==1)
     helpdlg('PRE-ICTAL');
     
     elseif (ani_security ==2)
     helpdlg('ICTAL ');

end 
end
%%






