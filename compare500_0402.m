% % % % %%
% % % % % % % % % % cd DeepCAD
% % % % % % cd N2N_WF500
cd n2n_SIM_ROIs
N2 = [];
files = dir('*.csv');
i = 1;
for file = files'
    csv = readtable(file.name);
    N2(i,:) = csv.Mean;
    i = i + 1;
end
cd ..

% cd pWF_ROIs_noDC
% WF = [];
% files = dir('*.csv');
% i = 1;
% for file = files'
%     csv = readtable(file.name);
%     WF(i,:) = csv.Mean;
%     i = i + 1;
% end
% cd ..
cd pWF_MC
WF = [];
files = dir('*.csv');
i = 1;
for file = files'
    csv = readtable(file.name);
    WF(i,:) = csv.Mean;
    i = i + 1;
end
cd ..

cd MC_pHiLo_DC
HN = [];
files = dir('*.csv');
i = 1;
for file = files'
    csv = readtable(file.name);
    HN(i,:) = csv.Mean;
    i = i + 1;
end
cd ..
% % 
cd MC_pHiLo_noDC
PH = [];
files = dir('*.csv');
i = 1;
for file = files'
    csv = readtable(file.name);
    PH(i,:) = csv.Mean;
    i = i + 1;
end
cd ..
% % 
% % 
% % 
% cd n2n_pWF_ROIs
% pN = [];
% files = dir('*.csv');
% i = 1;
% for file = files'
%     csv = readtable(file.name);
%     pN(i,:) = csv.Mean;
%     i = i + 1;
% end
% cd ..
cd pWF_N2N_MC
pN = [];
files = dir('*.csv');
i = 1;
for file = files'
    csv = readtable(file.name);
    pN(i,:) = csv.Mean;
    i = i + 1;
end
cd ..
% % 
% cd sim_noDC
% OS = [];
% files = dir('*.csv');
% i = 1;
% for file = files'
%     csv = readtable(file.name);
%     OS(i,:) = csv.Mean;
%     i = i + 1;
% end
% cd ..
cd SIM_SEED_MC_ROIs
OS = [];
files = dir('*.csv');
i = 1;
for file = files'
    csv = readtable(file.name);
    OS(i,:) = csv.Mean;
    i = i + 1;
end
cd ..







% % 
% 
N2 = N2(:,1:59988);
% % % 
% % % % 
% % % % % % 
% % % % % % 
% % % % % 
% % % % % 
% % % % % % 
% % % % % % % 
% % % % % % 

%%
%clf

noise = @(x) mean(x)/std(x);

% % %
ts = 1;
tf = 59988;
fr = 500;
xVals = linspace(0,((tf-ts)/fr)*1000,(tf-ts)+1);
xVals2 = linspace(0,((tf-ts)/fr)*1000,((tf-ts)+1)/10);

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



% Add noise to pWF
%pWF_WGN = zeros(size(WF));


p = 99;
for i = 1:6

    %pWF_WGN(i,:) = awgn(WF(i,:),43,'measured');

    % std_WF(1,i) = noise(WF(i,20000:3:22000));
    % std_N2(1,i) = noise(N2(i,20000:3:22000));
    % std_WG(1,i) = noise(pWF_WGN(i,20000:3:22000));

    dFF_SIM(i,:) = dFF2(N2(i,1:tf)', xVals);
    dFF_pWF(i,:) = dFF2(WF(i,1:tf)', xVals);
    dFF_pHL(i,:) = dFF2(PH(i,1:tf)', xVals);
    dFF_pHN(i,:) = dFF2(HN(i,1:tf)', xVals);
    dFF_OSS(i,:) = dFF2(OS(i,1:tf)', xVals);
    dFF_pWN(i,:) = dFF2(pN(i,1:tf)', xVals);


    %dFF_WG(i,:) = dFF2(pWF_WGN(i,:)', xVals);

    %peaks_pWF(1,i) = max(dFF_pWF(i,:));
    peaks_pWF(1,i) = prctile(dFF_pWF(i,:),p);
    std_pWF(1,i) = std(dFF_pWF(i,22250:24000));
    pnr_pWF(1,i) = peaks_pWF(1,i)/std_pWF(1,i);
    pnr_pWF_norm(1,i) = pnr_pWF(1,i)/pnr_pWF(1,i);

    %peaks_SIM(1,i) = max(dFF_SIM(i,:));
    peaks_SIM(1,i) = prctile(dFF_SIM(i,:),p);
    std_SIM(1,i) = std(dFF_SIM(i,22250:24000));
    pnr_SIM(1,i) = peaks_SIM(1,i)/std_SIM(1,i);
    pnr_SIM_norm(1,i) = pnr_SIM(1,i)/pnr_pWF(1,i);


    %peaks_pHL(1,i) = max(dFF_pHL(i,:));
    peaks_pHL(1,i) = prctile(dFF_pHL(i,:),p);
    std_pHL(1,i) = std(dFF_pHL(i,22250:24000));
    pnr_pHL(1,i) = peaks_pHL(1,i)/std_pHL(1,i);
    pnr_pHL_norm(1,i) = pnr_pHL(1,i)/pnr_pWF(1,i);

    %peaks_pHN(1,i) = max(dFF_pHN(i,:));
    peaks_pHN(1,i) = prctile(dFF_pHN(i,:),p);
    std_pHN(1,i) = std(dFF_pHN(i,22250:24000));
    pnr_pHN(1,i) = peaks_pHN(1,i)/std_pHN(1,i);
    pnr_pHN_norm(1,i) = pnr_pHN(1,i)/pnr_pWF(1,i);


    %peaks_OSS(1,i) = max(dFF_OSS(i,:));
    peaks_OSS(1,i) = prctile(dFF_OSS(i,:),p);
    std_OSS(1,i) = std(dFF_OSS(i,22250:24000));
    pnr_OSS(1,i) = peaks_OSS(1,i)/std_OSS(1,i);
    pnr_OSS_norm(1,i) = pnr_OSS(1,i)/pnr_pWF(1,i);

    %peaks_pWN(1,i) = max(dFF_pWN(i,:));
    peaks_pWN(1,i) = prctile(dFF_pWN(i,:),p);
    std_pWN(1,i) = std(dFF_pWN(i,22250:24000));
    pnr_pWN(1,i) = peaks_pWN(1,i)/std_pWN(1,i);
    pnr_pWN_norm(1,i) = pnr_pWN(1,i)/pnr_pWF(1,i);
    % 
    % 
    % 



    % peaks_WG(1,i) = max(dFF_WG(i,:));
    % std_WG(1,i) = std(dFF_WG(i,21000:22000));
    % pnr_WG(1,i) = peaks_WG(1,i)/std_WG(1,i);



    %HPF_OS(i,:) = highpass(dFF_SIM(i,:),15,100);


end
subplot(2,3,1)
bar(pnr_pWF_norm)
title("pWF sPNR (Mean = 21.78)")
xlabel("Cell #")
ylabel("sPNR")
PNR_AVG_pWF = mean(pnr_pWF_norm);

subplot(2,3,2)
bar(pnr_OSS_norm)
title("OS-SIM sPNR (Mean = 7.52)")
xlabel("Cell #")
ylabel("sPNR")
PNR_AVG_OSS = mean(pnr_OSS_norm);


subplot(2,3,3)
bar(pnr_pHL_norm)
PNR_AVG_pHL = mean(pnr_pHL_norm);
title("pHiLo sPNR (Mean = 28.69)")
xlabel("Cell #")
ylabel("sPNR")



subplot(2,3,4)
bar(pnr_pWN_norm)
PNR_AVG_pWN = mean(pnr_pWN_norm);
title("pWF-N2N sPNR (Mean = 23.42)")
xlabel("Cell #")
ylabel("sPNR")



subplot(2,3,5)
bar(pnr_SIM_norm)
PNR_AVG_SIM = mean(pnr_SIM_norm);
title("OS-SIM-N2N sPNR (Mean = 10.89)")
xlabel("Cell #")
ylabel("sPNR")

subplot(2,3,6)
bar(pnr_pHN_norm)
PNR_AVG_pHN = mean(pnr_pHN_norm);
title("pHiLo-N2N sPNR (Mean = 37.57)")
xlabel("Cell #")
ylabel("sPNR")
%%
clf
%x = ["pWF", "OS-SIM", "pHiLo"];
x = [1 2 3  4 5 6];
x2 = [1 2 3];
%data = [21.78, 7.52, 28.69, 23.42, 10.89, 37.57];
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

clr_noDC = [1 0 0; 1 0 0; 0 0 1;0 0 1;  0 1 0;0 1 0 ];
b2 = bar(x,mean_PNR_vals, 'facecolor', 'flat');     
b2.CData= clr_noDC;

% bar(x,mean_PNR_vals)                
yticks([1 2])
xticks([])
% hold on

hold on
er = errorbar(x,mean_PNR_vals,errs);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

% hold on
% er = errorbar(x2,mean_PNR_vals_noDC,errs_noDC);    
% er.Color = [0 0 0];                            
% er.LineStyle = 'none'; 
% 
% hold off

%%

mean_corr_vals = [mean_corr_pWF, mean_corr_pWN, mean_corr_OSS, mean_corr_SIM,mean_corr_pHL, mean_corr_pHN];

std_pWF_CC = std(pWF_Corrs);
std_OSS_CC = std(OSS_Corrs);
std_pHL_CC = std(pHL_Corrs);
std_pWN_CC = std(pWN_Corrs);
std_SIM_CC = std(SIM_Corrs);
std_pHN_CC = std(pHN_Corrs);


errs_CC = [std_pWF_CC, std_pWN_CC,  std_OSS_CC,  std_SIM_CC, std_pHL_CC, std_pHN_CC]/2;

b = bar(x,mean_corr_vals, FaceColor="flat");
b.CData = clr_noDC;

%yticks([0 15 30 45])
xticks([])
yticks([0.25 0.5])
hold on

er = errorbar(x,mean_corr_vals,errs_CC);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  


hold off






%%







%%



% % dFF_pWF = [dFF_pWF1,dFF_pWF2,dFF_pWF3,dFF_pWF4,dFF_pWF5,dFF_pWF6,dFF_pWF7,dFF_pWF8,dFF_pWF9,dFF_pWF10,dFF_pWF11];
% maxes_OS = zeros(1,11);
% maxes_pWF = zeros(1,11);
% % 
% % for ii = 1:11
% %     maxes_OS(1,ii) = max(dFF_SIM(:,ii));
% %     maxes_pWF(1,ii) = max(dFF_pWF(:,ii));
% % end
% % avg_max_SIM = mean(maxes_OS);
% % avg_max_pWF = mean(maxes_pWF);
% % dFF_ratio = avg_max_SIM/avg_max_pWF;
% % 
% % 
% % 
num_Cells = 6;
fr = 500;
offset_S = 1.3;
WF_scale = 8.0;
HL_scale = 1;
r_scale = 0.3;
xi = 22250*1000/fr; %Paper!
%xi = 33000*1000/fr;
%xf = 41000*1000/fr;
xf = 24000*1000/fr; %Paper!
%xf = 116000;
wi = 47000/10;
wf = 49000/10;
yi = 1;
yf = 7;
clf
for i = 1:5

    subplot(2,3,5)
    if i == 1
        plot(xVals,(smooth(dFF_SIM(i,:),3)+6*offset_S)/offset_S,LineWidth= 1.0)
    else
        plot(xVals,(smooth(dFF_SIM(i,:),3)+i*offset_S)/offset_S,LineWidth= 1.0)
    end 
    % plot(xVals, dFF_SIM(:,i)+i*offset_S-1,LineWidth= 1.0)
    xlim([xi,xf])
    ylim([yi,yf])
    %xticks([])
    yticks([])
    %ylabel("SIM")
    %xlabel("Time (ms)")
    ylabel("Cell #")
    hold on
    title("OS-SIM N2N (130% dF/F spacing)")

    subplot(2,3,1)
    if i == 1
        plot(xVals,(smooth(WF_scale*dFF_pWF(i,:),3)+6*offset_S)/offset_S,LineWidth= 1.0)
    else
        plot(xVals,(smooth(WF_scale*dFF_pWF(i,:),3)+i*offset_S)/offset_S,LineWidth= 1.0)
    end
    %plot(xVals, 2*dFF(dFF_n2n(10,ts:tf)',xVals)+i*offset_S)
    xlim([xi,xf])
    ylim([yi,yf])
    %xticks([])
    yticks([])
    ylabel("Cell #")
    xlabel("Time (ms)")
    hold on
    title("8.0 * pWF")

    subplot(2,3,3)
    if i == 1
        plot(xVals,(smooth(HL_scale*dFF_pHL(i,:),3)+6*offset_S)/offset_S,LineWidth= 1.0)
    else
        plot(xVals,(smooth(HL_scale*dFF_pHL(i,:),3)+i*offset_S)/offset_S,LineWidth= 1.0)
    end

    %plot(xVals, 2*dFF(dFF_n2n(10,ts:tf)',xVals)+i*offset_S)
    xlim([xi,xf])
    ylim([yi,yf])
    %xticks([])
    yticks([])
    ylabel("Cell #")
    xlabel("Time (ms)")
    hold on
    title("pHiLo")

    subplot(2,3,6)
    if i == 1
        plot(xVals,(smooth(HL_scale*dFF_pHN(i,:),3)+6*offset_S)/offset_S,LineWidth= 1.0)
    else
        plot(xVals,(smooth(HL_scale*dFF_pHN(i,:),3)+i*offset_S)/offset_S,LineWidth= 1.0)
    end
    %plot(xVals, 2*dFF(dFF_n2n(10,ts:tf)',xVals)+i*offset_S)
    xlim([xi,xf])
    ylim([yi,yf])
    %xticks([])
    yticks([])
    ylabel("Cell #")
    xlabel("Time (ms)")
    hold on
    title("pHiLo N2N")

    subplot(2,3,2)
    if i == 1
        plot(xVals,(smooth(1.0*dFF_OSS(i,:),3)+6*offset_S)/offset_S,LineWidth= 1.0)
    else
        plot(xVals,(smooth(1.0*dFF_OSS(i,:),3)+i*offset_S)/offset_S,LineWidth= 1.0)
    end
    %plot(xVals, 2*dFF(dFF_n2n(10,ts:tf)',xVals)+i*offset_S)
    xlim([xi,xf])
    ylim([yi,yf])
    %xticks([])
    yticks([])
    ylabel("Cell #")
    xlabel("Time (ms)")
    hold on
    title("OS-SIM")
    % 
    subplot(2,3,4)
    if i == 1
        plot(xVals,(smooth(WF_scale*dFF_pWN(i,:),3)+6*offset_S)/offset_S,LineWidth= 1.0)
    else
        plot(xVals,(smooth(WF_scale*dFF_pWN(i,:),3)+i*offset_S)/offset_S,LineWidth= 1.0)
    end

    %plot(xVals, 2*dFF(dFF_n2n(10,ts:tf)',xVals)+i*offset_S)
    xlim([xi,xf])
    ylim([yi,yf])
    %xticks([])
    yticks([])
    ylabel("Cell #")
    xlabel("Time (ms)")
    hold on
    title("8.0 * pWF N2N")


    % subplot(3,4,9)
    % plot(xVals,dFF_WG(i,:),LineWidth= 1.0)
    % 
    % %plot(xVals, 2*dFF(dFF_n2n(10,ts:tf)',xVals)+i*offset_S)
    % xlim([xi,xf])
    % %ylim([yi,yf])
    % ylabel("dF/F")
    % hold on
    % title("pWF + WGN dF/F")


end
% 
% % subplot(3,4,6)
% % histogram(dFF_SIM(1,:))
% % title("SIM dF/F Values")
% % set(gca,'yscale','log')
% % 
% % subplot(3,4,2)
% % histogram(dFF_pWF(1,:))
% % title("pWF dF/F Values")
% % set(gca,'yscale','log')
% % 
% % subplot(3,4,10)
% % histogram(dFF_WG(1,:))
% % title("pWF + WGN dF/F Values")
% % set(gca,'yscale','log')


%%
% % 
num_Cells = 6;

clf
corr_SIM = corrcoef(dFF_SIM(1:num_Cells,:)');
corr_pWF = corrcoef(dFF_pWF(1:num_Cells,:)');

%corr_WG = corrcoef(dFF_WG(1:num_Cells,:)');

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



%mean_corr_WG = mean(corr_WG, "all");
%ratio = mean_corr_OS/mean_corr_WF;   

SIM_Corrs = [];
pWF_Corrs = [];
OSS_Corrs = [];
pWN_Corrs = [];
pHL_Corrs = [];
pHN_Corrs = [];



%hm = zeros(num_Cells,num_Cells);
%hm2 = zeros(num_Cells,num_Cells);

for i = 1:num_Cells
    for l = 1:num_Cells
        if i < l 
            %hm(i,l) = corr_OS(i,l);
            %hm2(i,l) = corr_OS(i,l);
            SIM_Corrs(end+1) = corr_SIM(i,l);
            pWF_Corrs(end+1) = corr_pWF(i,l);
            OSS_Corrs(end+1) = corr_OSS(i,l);
            pWN_Corrs(end+1) = corr_pWN(i,l);
            pHL_Corrs(end+1) = corr_pHL(i,l);
            pHN_Corrs(end+1) = corr_pHN(i,l);





        elseif i>l
            %hm(i,l) = corr_WF(i,l);
            %hm2(i,l) = corr_WG(i,l);
            %wf_Corrs(end+1) = corr_WF(i,l);
        else 
            %hm(i,l) = NaN;
            %hm2(i,l) = NaN;
        end
    end
end

subplot(2,3,1)
histogram(pWF_Corrs)
title("pWF Corr. Coefs. (Mean = 0.523)")
xlabel("Correlation Coeff.")
ylabel("Occurances")

subplot(2,3,2)
histogram(OSS_Corrs)
title("OS-SIM Corr. Coefs. (Mean = 0.179)")
xlabel("Correlation Coeff.")
ylabel("Occurances")

subplot(2,3,3)
histogram(pHL_Corrs)
title("pHiLo Corr. Coefs. (Mean = 0.199)")
xlabel("Correlation Coeff.")
ylabel("Occurances")

subplot(2,3,4)
histogram(pWN_Corrs)
title("pWF-N2N Corr. Coefs. (Mean = 0.527)")
xlabel("Correlation Coeff.")
ylabel("Occurances")

subplot(2,3,5)
histogram(SIM_Corrs)
title("OS-SIM-N2N Corr. Coefs. (Mean = 0.180)")
xlabel("Correlation Coeff.")
ylabel("Occurances")

subplot(2,3,6)
histogram(pHN_Corrs)
title("pHiLo-N2N Corr. Coefs. (Mean = 0.198)")
xlabel("Correlation Coeff.")
ylabel("Occurances")

% 
% 
% 
% 
% 






% subplot(2,4,3)
% % % %heatmap(corr_OS,'ColorLimits',[0.1 0.85], 'Colormap', jet)
% % % %heatmap(corr_OS,'ColorLimits',[0.0 1])
% heatmap(hm,'ColorLimits',[0.0 1])
% ylabel("Cell # pWF")
% xlabel("Cell # pWF")
% title("Heatmap with pWF")
% 
% subplot(2,4,7)
% % % %heatmap(corr_OS,'ColorLimits',[0.1 0.85], 'Colormap', jet)
% % % %heatmap(corr_OS,'ColorLimits',[0.0 1])
% heatmap(hm2,'ColorLimits',[0.0 1])
% ylabel("Cell # pWF+WGN")
% xlabel("Cell # pWF+WGN")
% title("Heatmap with pWF+WGN")
% 
% 
% subplot(2,4,4)
% y = [mean_corr_OS, mean_corr_WF, mean_corr_WG];
% x = ["OS" "WF" "WF+WG"];
% bar(x,y)
% title("Mean Correlation Coefficient")
% 
% subplot(2,4,8)
% y = [PNR_AVG_OS, PNR_AVG_WF, PNR_AVG_WG ];
% x = ["OS" "WF" "WF+WG"];
% bar(x,y)
% title("Mean Spike Peak-to-Noise Ratio")


%plot(xVals, dFF_SIM(:,10))
%xlim([46000,48000])
% % % % 
% % % 
% % % 
% % % 
% % % % 
% subplot(3,3,3)
% hist(wf_Corrs)
% xlabel("Corr. Coeff. (N2N)")
% %xlim([0,0.35])
% title("Noise2Noise AVG. Coeff = 0.2357")
% 
% subplot(3,3,6)
% hist(os_Corrs)
% xlabel("Corr. Coeff. (FB)")
% %xlim([0,0.35])
% title("Fish Brains AVG. Coeff = 0.2668")


%%
clf
% j = 0;
% k = 0;
% l = 0;
% jj = 0;
% kk = 0;
% ll = 0;
% 
% hits_SIM = zeros(tf,1);
% hits_pWN = zeros(tf,1);
% hits_pHN = zeros(tf,1);
% 
% hits_OSS = zeros(tf,1);
% hits_pWF = zeros(tf,1);
% hits_pHL = zeros(tf,1);
% 
% thresh_scale = 0.3;
% 
% % 
% for z = 15000:59994
%     for i = 1:6
%         if dFF_pWF(i,z) >= max(dFF_pWF(i,:))*thresh_scale
%             j = j + 1;
%         end
%         hits_pWF(z,1) = j;
% 
%         if dFF_OSS(i,z) >= max(dFF_OSS(i,:))*thresh_scale
%             k = k + 1;
%         end
%         hits_OSS(z,1) = k;
% 
%         if dFF_pHL(i,z) >= max(dFF_pHL(i,:))*thresh_scale
%             l = l + 1;
%         end
%         hits_pHL(z,1) = l;
% 
% 
% 
%         if dFF_pWN(i,z) >= max(dFF_pWN(i,:))*thresh_scale
%             jj = jj + 1;
%         end
%         hits_pWN(z,1) = jj;
% 
%         if dFF_SIM(i,z) >= max(dFF_SIM(i,:))*thresh_scale
%             kk = kk + 1;
%         end
%         hits_SIM(z,1) = kk;
% 
%         if dFF_pHN(i,z) >= max(dFF_pHN(i,:))*thresh_scale
%             ll = ll + 1;
%         end
%         hits_pHN(z,1) = ll;
% 
% 
% 
%     end
% end
% % 
% % 
% % 
% % 
% 



%subplot(1,3,3)
%plot(hits_pWF(:,1),'r',LineWidth=2)
%hold on 
% plot(hits_OSS(:,1),'b',LineWidth=2)
% plot(hits_pHL(:,1),'g',LineWidth=2)

xVals_hits = xVals - 30000;

plot(xVals_hits, hits_pWN(:,1),'r',LineWidth=2)
hold on 
plot(xVals_hits, hits_SIM(:,1),'b',LineWidth=2)
plot(xVals_hits, hits_pHN(:,1),'g',LineWidth=2)




xlabel("Time (ms)")
legend("pWF",'OS-SIM','pHiLo', Location = 'northwest')
xlim([xi-30000,xf-30000])
%ylim([yi,yf])
%xticks([])
%yticks([])
%ylim([0, max(hits_pWF(:,1))])
percent = 100*k/j;
ylabel('Frames above 30% Threshold')
% 
% percent_OS = 100*k/(8*29501);
% percent_WF = 100*j/(8*29501);%%

%%

% % t = Tiff('n2n_pHiLo.tif', 'r');
% % imageData = read_file(t);
% 
% %filename_hil = 'pHiLo_0point1_n2n.tif';
% %[hilData, ~, ~, ~] = read_file(filename_hil);
% 
% 
% 
% 
% 
% filename_SIM = 'NEIL_1_n90_500_b4_MMStack.ome_interleaved_NOmatchSIM_Reconstruction.tif';
% %filename_pWF = 'NEIL_1_p1_E_05_Iter_6254_outp_i_interleaved_NOmatchpWF_Reconstruction.tif';
% %filename_hil = 'n2n_pHiLo_16b.tif';
% 
% % t = Tiff('n2n_pHiLo.tif', 'r');
% % imageData = read(t);
% 
% 
% info = imfinfo(filename_pWF);
% numframe = 59994;
% for K = 15000 : numframe
%     rawframes_SIM(:,:,K-14999) = imread(filename_SIM, K);
%     %rawframes_pWF(:,:,K-14999) = imread(filename_pWF, K);
%     %rawframes_hil(:,:,K-14999) = imread(filename_hil, K);
% 
%     %rawframes_rSIM(:,:,K) = imread(filename_rSIM, K);
% 
% end
% simframes = mat2gray(rawframes_SIM);
% 
% 
% %pWFframes = mat2gray(rawframes_pWF);
% %rsimframes = mat2gray(rawframes_rSIM);
% 
% %hilframes = mat2gray(hilData(:,:,15000:59994));
% 
% % 
% corr_image_SIM = zeros(106,166);
% %corr_image_pWF = zeros(106,166);
% %corr_image_rSIM = zeros(84,230);
% %corr_image_hil = zeros(106,166);
% 
% full_os = fullOS.Mean;
% f_os = full_os(15000:59994);
% 
% % full_wf = FullpWF.Mean;
% % f_wf = full_wf(15000:59994);
% % 
% % full_hl = FullHL.Mean;
% % f_hl = full_hl(15000:59994);
% 
% 
% for x=1:166
%     for y = 1:106
%         cc_sim = corrcoef(simframes(y,x,:),f_os);
%         %cc_pWF = corrcoef(pWFframes(y,x,:),f_wf);
%         %cc_rsim = corrcoef(rsimframes(y,x,:),rsim(1,:));
%         %cc_hil = corrcoef(hilframes(y,x,:),f_hl);
% 
%         corr_image_sim(y,x) = cc_sim(2,1);
%         %corr_image_pwf(y,x) = cc_pWF(2,1);
%         %corr_image_rsim(y,x) = cc_rsim(2,1);
%         %corr_image_hil(y,x) = cc_hil(2,1);
% 
% 
%     end
% end

clims = [0 1];
subplot(1,3,1)
imagesc(corr_image_pwf, clims)
colorbar
title("pWF Correlation Map")
xticks([])
yticks([])

subplot(1,3,2)
imagesc(corr_image_sim, clims)
colorbar
xticks([])
yticks([])
title("OS-SIM Correlation Map")

subplot(1,3,3)
imagesc(corr_image_hil, clims)
colorbar
xticks([])
yticks([])
title("HiLo Correlation Map")

% subplot(4,5,15)
% imagesc(corr_image_rsim)
% colorbar
% title("rSIM Correlation Map")
% 
% subplot(4,5,20)
% imagesc(corr_image_hil)
% colorbar
% title("pHilo Correlation Map")

%%

%%

% clf
% % y = sim10.Mean(ts:tf,1);
% opol = 15;
% [p,s,mu] = polyfit(xVals',y(:,1),opol);
% f_y = polyval(p,xVals',[],mu);
% 
% dt_ecgnl = y - f_y;
% %dFF = dt_ecgnl./f_y;

perc(1.41,1.84)

function percent_increase = perc(x,denoised_x)

    percent_increase = (denoised_x-x)/x;

end



function [dF_F,F0] = dFF(y,xVals)

    % y2 = msbackadj(xVals', ...              % transpose to column vector
    %                y, ...              % transpose to column vector
    %                'WindowSize', 1000, ...
    %                'StepSize', 1000, ...
    %                'Quantile', 0.1, ...
    %                'ShowPlot',  'no');
    % 
    % F0 = y'-y2';
    % dF_F = y2./F0';
    dF_F = y./prctile(y(:,1),[10],"all");
end

function [dF_F,F0] = dFF2(y,xVals)

    y2 = msbackadj(xVals', ...              % transpose to column vector
                   y, ...              % transpose to column vector
                   'WindowSize', 5000, ...
                   'StepSize', 5000, ...
                   'Quantile', 0.1, ...
                   'ShowPlot',  'no');

    F0 = y'-y2';
    dF_F = y2./F0';
   
end

function [imgstack, xPixels, yPixels, tPoints] = read_file(fileName)

    Data = bfopen(fileName);

    val = Data{1,1};

    val2 = val(:,1);

    [xPixels,yPixels] = size(val2{1,1});

    tPoints = length(val2);

    imgstack = zeros(xPixels, yPixels, tPoints, 'uint16');

    for i=1:tPoints

        imgstack(:,:,i)=cell2mat(val2(i));

    end

   

end

% function snr = SNR(QE,S,Ib,Nr)
% 
%     snr = (QE*S)
% 
% end 








