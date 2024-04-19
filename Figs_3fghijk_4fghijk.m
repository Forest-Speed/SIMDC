% Import time courses
% 
% 
% cd analyze_100
% 
% cd ROIs_N2N
% N2 = [];
% files = dir('*.csv');
% i = 1;
% for file = files'
%     csv = readtable(file.name);
%     N2(i,:) = csv.Mean;
%     i = i + 1;
% end
% cd ..
% 
% cd ROIs_pWF_noDEEP
% WF = [];
% files = dir('*.csv');
% i = 1;
% for file = files'
%     csv = readtable(file.name);
%     WF(i,:) = csv.Mean;
%     i = i + 1;
% end
% cd ..
% 
% cd ROIs_pHiLo
% HN = [];
% files = dir('*.csv');
% i = 1;
% for file = files'
%     csv = readtable(file.name);
%     HN(i,:) = csv.Mean;
%     i = i + 1;
% end
% cd ..
% 
% cd ROIs_pHiLo_noDEEP
% PH = [];
% files = dir('*.csv');
% i = 1;
% for file = files'
%     csv = readtable(file.name);
%     PH(i,:) = csv.Mean;
%     i = i + 1;
% end
% cd ..
% 
% 
% 
% cd ROIs_pWF_N2N
% pN = [];
% files = dir('*.csv');
% i = 1;
% for file = files'
%     csv = readtable(file.name);
%     pN(i,:) = csv.Mean;
%     i = i + 1;
% end
% cd ..
% 
% cd ROIs_OS_noDeep
% OS = [];
% files = dir('*.csv');
% i = 1;
% for file = files'
%     csv = readtable(file.name);
%     OS(i,:) = csv.Mean;
%     i = i + 1;
% end
% cd ..
% 
% 
% OS = OS(:,1:35994);
%%
dFF_SIM = zeros(size(N2));
dFF_pWF = zeros(size(N2));
dFF_WG = zeros(size(N2));
HPF_OS = zeros(size(N2));

dFF_pHL = zeros(size(N2));
dFF_pHN = zeros(size(N2));

dFF_OSS = zeros(size(N2));
dFF_pWN = zeros(size(N2));



std_WF = zeros(1,6);
std_N2 = zeros(1,6);
std_WG = zeros(1,6);
peaks_OS = zeros(1,6);
std_OS = zeros(1,6);
pnr_OS = zeros(1,6);
peaks_WG = zeros(1,6);
pnr_WG = zeros(1,6);
peaks_WF = zeros(1,6);
pnr_WF = zeros(1,6);

pnr_pWF_norm = zeros(1,6);
pnr_OSS_norm = zeros(1,6);
pnr_pHL_norm = zeros(1,6);
pnr_pWN_norm = zeros(1,6);
pnr_SIM_norm = zeros(1,6);
pnr_pHN_norm = zeros(1,6);



n_ts = 210000;
n_tf = 230000;
ns = n_ts/100;
nf = n_tf/100;


p = 99;
for i = 1:6

    dFF_SIM(i,:) = dFF2(N2(i,:)', xVals);
    dFF_pWF(i,:) = dFF2(WF(i,:)', xVals);
    dFF_pHL(i,:) = dFF2(PH(i,:)', xVals);
    dFF_pHN(i,:) = dFF2(HN(i,:)', xVals);
    dFF_OSS(i,:) = dFF2(OS(i,:)', xVals);
    dFF_pWN(i,:) = dFF2(pN(i,:)', xVals);




    %peaks_pWF(1,i) = max(dFF_pWF(i,:));
    peaks_pWF(1,i) = prctile(dFF_pWF(i,:),p);
    std_pWF(1,i) = std(dFF_pWF(i,ns:nf));
    pnr_pWF(1,i) = peaks_pWF(1,i)/std_pWF(1,i);
    pnr_pWF_norm(1,i) = pnr_pWF(1,i)/pnr_pWF(1,i);

    %peaks_SIM(1,i) = max(dFF_SIM(i,:));
    peaks_SIM(1,i) = prctile(dFF_SIM(i,:),p);
    std_SIM(1,i) = std(dFF_SIM(i,ns:nf));
    pnr_SIM(1,i) = peaks_SIM(1,i)/std_SIM(1,i);
    pnr_SIM_norm(1,i) = pnr_SIM(1,i)/pnr_pWF(1,i);

    %peaks_pHL(1,i) = max(dFF_pHL(i,:));
    peaks_pHL(1,i) = prctile(dFF_pHL(i,:),p);
    std_pHL(1,i) = std(dFF_pHL(i,ns:nf));
    pnr_pHL(1,i) = peaks_pHL(1,i)/std_pHL(1,i);
    pnr_pHL_norm(1,i) = pnr_pHL(1,i)/pnr_pWF(1,i);

    %peaks_pHN(1,i) = max(dFF_pHN(i,:));
    peaks_pHN(1,i) = prctile(dFF_pHN(i,:),p);
    std_pHN(1,i) = std(dFF_pHN(i,ns:nf));
    pnr_pHN(1,i) = peaks_pHN(1,i)/std_pHN(1,i);
    pnr_pHN_norm(1,i) = pnr_pHN(1,i)/pnr_pWF(1,i);


    %peaks_OSS(1,i) = max(dFF_OSS(i,:));
    peaks_OSS(1,i) = prctile(dFF_OSS(i,:),p);
    std_OSS(1,i) = std(dFF_OSS(i,ns:nf));
    pnr_OSS(1,i) = peaks_OSS(1,i)/std_OSS(1,i);
    pnr_OSS_norm(1,i) = pnr_OSS(1,i)/pnr_pWF(1,i);

    %peaks_pWN(1,i) = max(dFF_pWN(i,:));
    peaks_pWN(1,i) = prctile(dFF_pWN(i,:),p);
    std_pWN(1,i) = std(dFF_pWN(i,ns:nf));
    pnr_pWN(1,i) = peaks_pWN(1,i)/std_pWN(1,i);
    pnr_pWN_norm(1,i) = pnr_pWN(1,i)/pnr_pWF(1,i);




end






%% Fig 3f-h and 4f-h

ts = 1;
tf = 35994;
fr = 100;
xVals = linspace(0,((tf-ts)/fr)*1000,(tf-ts)+1);

num_Cells = 6;
fr = 100;
offset_S = 1.5;
WF_scale = 5.0;
HL_scale = 0.55;
r_scale = 0.3;
xi = 120000; 
xf = 360000;



wi = 47000/10;
wf = 49000/10;
yi = 0.5;
yf = 6.9;
clf
for i = 1:6

    subplot(2,3,5)
    if i <= 10
        plot(xVals,(smooth(dFF_SIM(i,:),3)+i*offset_S)/offset_S,LineWidth= 1.0)
    else
        x3 = linspace(wi,wf,100);
        y3 = ones(100,1)*10;
        plot(x3, y3,'Color','magenta',LineWidth= 5.0)
    end 
    xlim([xi,xf])
    ylim([yi,yf])
    ylabel("Cell #")
    hold on
    title("OS-SIM N2N")

    subplot(2,3,1)
    plot(xVals,(smooth(WF_scale*dFF_pWF(i,:),3)+i*offset_S)/offset_S,LineWidth= 1.0)
    xlim([xi,xf])
    ylim([yi,yf])
    ylabel("Cell #")
    xlabel("Time (ms)")
    hold on
    title("5.0 * pWF")

    subplot(2,3,3)
    plot(xVals,(smooth(HL_scale*dFF_pHL(i,:),3)+i*offset_S)/offset_S,LineWidth= 1.0)
    xlim([xi,xf])
    ylim([yi,yf])
    ylabel("Cell #")
    xlabel("Time (ms)")
    hold on
    title("0.55 * pHiLo")

    subplot(2,3,6)
    plot(xVals,(smooth(HL_scale*dFF_pHN(i,:),3)+i*offset_S)/offset_S,LineWidth= 1.0)
    xlim([xi,xf])
    ylim([yi,yf])
    ylabel("Cell #")
    xlabel("Time (ms)")
    hold on
    title("0.55 * pHiLo N2N")

    subplot(2,3,2)
    plot(xVals,(smooth(2.0*dFF_OSS(i,:),3)+i*offset_S)/offset_S,LineWidth= 1.0)
    xlim([xi,xf])
    ylim([yi,yf])
    ylabel("Cell #")
    xlabel("Time (ms)")
    hold on
    title("2.0 * OS-SIM")

    subplot(2,3,4)
    plot(xVals,(smooth(2.5*WF_scale*dFF_pWN(i,:),3)+i*offset_S)/offset_S,LineWidth= 1.0)
    xlim([xi,xf])
    ylim([yi,yf])
    ylabel("Cell #")
    xlabel("Time (ms)")
    hold on
    title("7.5 * pWF N2N")


end


%% Figs 3I and 4I

clf 

PNR_AVG_pWF = mean(pnr_pWF_norm);
PNR_AVG_OSS = mean(pnr_OSS_norm);
PNR_AVG_pHL = mean(pnr_pHL_norm);
PNR_AVG_pWN = mean(pnr_pWN_norm);
PNR_AVG_SIM = mean(pnr_SIM_norm);
PNR_AVG_pHN = mean(pnr_pHN_norm);


x = [1 2 3  4 5 6];
x2 = [1 2 3];

mean_PNR_vals = [PNR_AVG_pWF, PNR_AVG_pWN, PNR_AVG_OSS, PNR_AVG_SIM, PNR_AVG_pHL, PNR_AVG_pHN];
mean_PNR_vals_noDC = [PNR_AVG_pWF, PNR_AVG_OSS,  PNR_AVG_pHL];
mean_PNR_vals_DC = [PNR_AVG_pWN, PNR_AVG_SIM,  PNR_AVG_pHN];


std_pWF_PNR = std(pnr_pWF_norm);
std_OSS_PNR = std(pnr_OSS_norm);
std_pHL_PNR = std(pnr_pHL_norm);
std_pWN_PNR = std(pnr_pWN_norm);
std_SIM_PNR = std(pnr_SIM_norm);
std_pHN_PNR = std(pnr_pHN_norm);


errs = [std_pWF_PNR, std_pWN_PNR,  std_OSS_PNR,  std_SIM_PNR, std_pHL_PNR, std_pHN_PNR]/2;
errs_noDC = [std_pWF_PNR,   std_OSS_PNR,  std_pHL_PNR]/2;
errs_DC = [std_pWN_PNR, std_SIM_PNR, std_pHN_PNR]/2;


clr = [1 0 0; 1 0 0; 0 0 1; 0 0 1; 0 1 0; 0 1 0];
clr_noDC = [1 0 0; 0 0 1;  0 1 0 ];


subplot(1,2,1)
b1 = bar(x2,mean_PNR_vals_noDC, 'facecolor', 'flat');     
b1.CData= clr_noDC;
hold on
er = errorbar(x2,mean_PNR_vals_noDC,errs_noDC);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  



subplot(1,2,2)
b2 = bar(x2,mean_PNR_vals_DC, 'facecolor', 'flat');     
b2.CData= clr_noDC;
yticks([1 2 3])
xticks([])
hold on
er = errorbar(x2,mean_PNR_vals_DC,errs_DC);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  



%% Figs 3j and 4j

corr_SIM = corrcoef(dFF_SIM(1:num_Cells,:)');
corr_pWF = corrcoef(dFF_pWF(1:num_Cells,:)');
corr_OSS = corrcoef(dFF_OSS(1:num_Cells,:)');
corr_pWN = corrcoef(dFF_pWN(1:num_Cells,:)');
corr_pHL = corrcoef(dFF_pHL(1:num_Cells,:)');
corr_pHN = corrcoef(dFF_pHN(1:num_Cells,:)');


mean_corr_SIM = mean(corr_SIM, "all"); 
mean_corr_pWF = mean(corr_pWF, "all"); 
mean_corr_OSS = mean(corr_OSS, "all"); 
mean_corr_pWN = mean(corr_pWN, "all"); 
mean_corr_pHL = mean(corr_pHL, "all"); 
mean_corr_pHN = mean(corr_pHN, "all"); 


SIM_Corrs = [];
pWF_Corrs = [];
OSS_Corrs = [];
pWN_Corrs = [];
pHL_Corrs = [];
pHN_Corrs = [];

for i = 1:num_Cells
    for l = 1:num_Cells
        if i < l 
            SIM_Corrs(end+1) = corr_SIM(i,l);
            pWF_Corrs(end+1) = corr_pWF(i,l);
            OSS_Corrs(end+1) = corr_OSS(i,l);
            pWN_Corrs(end+1) = corr_pWN(i,l);
            pHL_Corrs(end+1) = corr_pHL(i,l);
            pHN_Corrs(end+1) = corr_pHN(i,l);
        end
    end
end

mean_corr_vals_noDC = [mean_corr_pWF,  mean_corr_OSS, mean_corr_pHL];
mean_corr_vals_DC = [mean_corr_pWN, mean_corr_SIM, mean_corr_pHN];

std_pWF_CC = std(pWF_Corrs);
std_OSS_CC = std(OSS_Corrs);
std_pHL_CC = std(pHL_Corrs);
std_pWN_CC = std(pWN_Corrs);
std_SIM_CC = std(SIM_Corrs);
std_pHN_CC = std(pHN_Corrs);


errs_CC = [std_pWF_CC, std_pWN_CC,  std_OSS_CC,  std_SIM_CC, std_pHL_CC, std_pHN_CC]/2;
errs_CC_noDC = [std_pWF_CC,  std_OSS_CC,   std_pHL_CC]/2;
errs_CC_DC = [ std_pWN_CC,    std_SIM_CC, std_pHN_CC]/2;

subplot(1,2,1)
b = bar(x2,mean_corr_vals_noDC, 'facecolor', 'flat');
b.CData= clr_noDC;
yticks([0 0.3 0.6])
xticks([])
hold on
er = errorbar(x2,mean_corr_vals_noDC,errs_CC_noDC);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 


subplot(1,2,2)
b = bar(x2,mean_corr_vals_DC, 'facecolor', 'flat');
b.CData= clr_noDC;
yticks([0 0.3 0.6])
xticks([])
hold on
er = errorbar(x2,mean_corr_vals_DC,errs_CC_DC);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

%% Figs 3k and 4k

clf
j = 0;
k = 0;
l = 0;
jj = 0;
kk = 0;
ll = 0;

hits_SIM = zeros(tf,1);
hits_pWN = zeros(tf,1);
hits_pHN = zeros(tf,1);

hits_OSS = zeros(tf,1);
hits_pWF = zeros(tf,1);
hits_pHL = zeros(tf,1);

thresh_scale = 0.3;


for z = 12000:35994
    for i = 1:6
        if dFF_pWF(i,z) >= max(dFF_pWF(i,12000:35994))*thresh_scale
            j = j + 1;
        end
        hits_pWF(z,1) = j;

        if dFF_OSS(i,z) >= max(dFF_OSS(i,12000:35994))*thresh_scale
            k = k + 1;
        end
        hits_OSS(z,1) = k;

        if dFF_pHL(i,z) >= max(dFF_pHL(i,12000:35994))*thresh_scale
            l = l + 1;
        end
        hits_pHL(z,1) = l;



        if dFF_pWN(i,z) >= max(dFF_pWN(i,12000:35994))*thresh_scale
            jj = jj + 1;
        end
        hits_pWN(z,1) = jj;

        if dFF_SIM(i,z) >= max(dFF_SIM(i,12000:35994))*thresh_scale
            kk = kk + 1;
        end
        hits_SIM(z,1) = kk;

        if dFF_pHN(i,z) >= max(dFF_pHN(i,12000:35994))*thresh_scale
            ll = ll + 1;
        end
        hits_pHN(z,1) = ll;



    end
end

frames = linspace(1,35994-12000,35995-12000);



subplot(1,2,1)
plot(frames, hits_pWF(12000:35994,1),'r',LineWidth=2)
hold on 
plot(frames,  hits_OSS(12000:35994,1),'b',LineWidth=2)
plot(frames, hits_pHL(12000:35994,1),'g',LineWidth=2)

subplot(1,2,2)
plot(frames, hits_pWN(12000:35994,1),'r',LineWidth=2)
hold on 
plot(frames, hits_SIM(12000:35994,1),'b',LineWidth=2)
plot(frames, hits_pHN(12000:35994,1),'g',LineWidth=2)




xlabel("Frames")
legend("pWF",'OS-SIM','pHiLo', Location = 'northwest')
percent = 100*k/j;
ylabel('Frames above 30 % Threshold')













%%
function [dF_F,F0] = dFF2(y,xVals)

    y2 = msbackadj(xVals', ...              % transpose to column vector
                   y, ...              % transpose to column vector
                   'WindowSize', 9000, ...
                   'StepSize', 9000, ...
                   'Quantile', 0.1, ...
                   'ShowPlot',  'no');

    F0 = y'-y2';
    dF_F = y2./F0';
   
end











