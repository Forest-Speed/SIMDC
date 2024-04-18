% % % % %
% % % % % % % % % cd DeepCAD
% %cd N2N_WF500
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
% 
% % 
% % % 
% % 
% 
% 
% % 
% % % 
% % 

% %%
% clf
% 
% noise = @(x) mean(x)/std(x);
% 
% % % %
ts = 1;
tf = 35994;
fr = 100;
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

pnr_pWF_norm = zeros(1,6);
pnr_OSS_norm = zeros(1,6);
pnr_pHL_norm = zeros(1,6);
pnr_pWN_norm = zeros(1,6);
pnr_SIM_norm = zeros(1,6);
pnr_pHN_norm = zeros(1,6);

% Add noise to pWF
pWF_WGN = zeros(size(WF));


n_ts = 210000;
n_tf = 230000;
ns = n_ts/100;
nf = n_tf/100;


p = 99;
for i = 1:6

    %pWF_WGN(i,:) = awgn(WF(i,:),43,'measured');

    % std_WF(1,i) = noise(WF(i,20000:3:22000));
    % std_N2(1,i) = noise(N2(i,20000:3:22000));
    % std_WG(1,i) = noise(pWF_WGN(i,20000:3:22000));

    dFF_SIM(i,:) = dFF2(N2(i,:)', xVals);
    dFF_pWF(i,:) = dFF2(WF(i,:)', xVals);
    dFF_pHL(i,:) = dFF2(PH(i,:)', xVals);
    dFF_pHN(i,:) = dFF2(HN(i,:)', xVals);
    dFF_OSS(i,:) = dFF2(OS(i,:)', xVals);
    dFF_pWN(i,:) = dFF2(pN(i,:)', xVals);


    %dFF_WG(i,:) = dFF2(pWF_WGN(i,:)', xVals);

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






    % peaks_WG(1,i) = max(dFF_WG(i,:));
    % std_WG(1,i) = std(dFF_WG(i,21000:22000));
    % pnr_WG(1,i) = peaks_WG(1,i)/std_WG(1,i);



    %HPF_OS(i,:) = highpass(dFF_SIM(i,:),15,100);


end
% subplot(2,3,1)
% %bar(pnr_pWF)
% bar(pnr_pWF_norm)
% PNR_AVG_pWF = mean(pnr_pWF_norm);
% title("pWF sPNR (Mean = " + PNR_AVG_pWF + ")")
% xlabel("Cell #")
% ylabel("sPNR")
% 
% 
% subplot(2,3,2)
% %bar(pnr_OSS)
% bar(pnr_OSS_norm)
% PNR_AVG_OSS = mean(pnr_OSS_norm);
% title("OS-SIM sPNR (Mean = " + PNR_AVG_OSS + ")")
% xlabel("Cell #")
% ylabel("sPNR")
% 
% 
% 
% subplot(2,3,3)
% %bar(pnr_pHL)
% bar(pnr_pHL_norm)
% PNR_AVG_pHL = mean(pnr_pHL_norm);
% title("pHiLo sPNR (Mean = " + PNR_AVG_pHL + ")")
% xlabel("Cell #")
% ylabel("sPNR")
% 
% 
% 
% subplot(2,3,4)
% %bar(pnr_pWN)
% bar(pnr_pWN_norm)
% PNR_AVG_pWN = mean(pnr_pWN_norm);
% title("pWF-DN sPNR (Mean = " + PNR_AVG_pWN + ")")
% xlabel("Cell #")
% ylabel("sPNR")
% 
% 
% 
% subplot(2,3,5)
% %bar(pnr_SIM)
% bar(pnr_SIM_norm)
% PNR_AVG_SIM = mean(pnr_SIM_norm);
% title("OS-SIM-DN sPNR (Mean = " + PNR_AVG_SIM + ")")
% xlabel("Cell #")
% ylabel("sPNR")
% 
% subplot(2,3,6)
% %bar(pnr_pHN)
% bar(pnr_pHN_norm)
% PNR_AVG_pHN = mean(pnr_pHN_norm);
% title("pHiLo-DN sPNR (Mean = " + PNR_AVG_pHN + ")")
% xlabel("Cell #")
% ylabel("sPNR")
% 
% 
% 
% % % WF_SNR_AVG = mean(std_WF);
% % PNR_AVG_OS = mean(pnr_OS);
% % PNR_AVG_WG = mean(pnr_WG);
% % PNR_AVG_WF = mean(pnr_WF);
% % 
% % subplot(3,3,1)
% % bar(pnr_WF)
% % 
% % subplot(3,3,4)
% % bar(pnr_OS)
% 
% % subplot(3,3,7)
% % bar(pnr_WG)
% 
% %% Magnitude Amplification
% 
% mag_pWF = mean(peaks_pWF);
% mag_OSS = mean(peaks_OSS);
% mag_pHL = mean(peaks_pHL);
% mag_pWN = mean(peaks_pWN);
% mag_SIM = mean(peaks_SIM);
% mag_pHN = mean(peaks_pHN);
% 
% mean_mag_noDC = [mag_pWF mag_OSS mag_pHL]/mag_pWF;
% mean_mag_DC = [mag_pWN mag_SIM mag_pHN]/mag_pWN;
% 
% 
% b2 = bar(x2,mean_mag_DC, 'facecolor', 'flat');     
% b2.CData= clr_noDC;
% 
% yticks([3 6])
% xticks([])
% 
% 
% 
% 
% 
% 
% 
% 
% %%
% clf
% %x = ["pWF", "OS-SIM", "pHiLo"];
% x = [1 2 3  4 5 6];
% x2 = [1 2 3];
% %data = [21.78, 7.52, 28.69, 23.42, 10.89, 37.57];
% mean_PNR_vals = [PNR_AVG_pWF, PNR_AVG_pWN, PNR_AVG_OSS, PNR_AVG_SIM, PNR_AVG_pHL, PNR_AVG_pHN];
% 
% mean_PNR_vals_noDC = [PNR_AVG_pWF, PNR_AVG_OSS,  PNR_AVG_pHL];
% mean_PNR_vals_DC = [PNR_AVG_pWN, PNR_AVG_SIM,  PNR_AVG_pHN];
% 
% 
% std_pWF_PNR = std(pnr_pWF_norm);
% std_OSS_PNR = std(pnr_OSS_norm);
% std_pHL_PNR = std(pnr_pHL_norm);
% std_pWN_PNR = std(pnr_pWN_norm);
% std_SIM_PNR = std(pnr_SIM_norm);
% std_pHN_PNR = std(pnr_pHN_norm);
% 
% 
% errs = [std_pWF_PNR, std_pWN_PNR,  std_OSS_PNR,  std_SIM_PNR, std_pHL_PNR, std_pHN_PNR]/2;
% errs_noDC = [std_pWF_PNR,   std_OSS_PNR,  std_pHL_PNR]/2;
% errs_DC = [std_pWN_PNR, std_SIM_PNR, std_pHN_PNR]/2;
% 
% 
% clr = [1 0 0; 1 0 0; 0 0 1; 0 0 1; 0 1 0; 0 1 0];
% clr_noDC = [1 0 0; 0 0 1;  0 1 0 ];
% 
% % b = bar(x,mean_PNR_vals, 'facecolor', 'flat');     
% % b.CData= clr;
% 
% % b1 = bar(x2,mean_PNR_vals_noDC, 'facecolor', 'flat');     
% % b1.CData= clr_noDC;
% 
% b2 = bar(x2,mean_PNR_vals_DC, 'facecolor', 'flat');     
% b2.CData= clr_noDC;
% 
% 
% 
% yticks([1 2 3])
% xticks([])
% hold on
% 
% er = errorbar(x2,mean_PNR_vals_DC,errs_DC);    
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% 
% 
% hold off
% 
% %%
% 
% mean_corr_vals = [mean_corr_pWF, mean_corr_pWN, mean_corr_OSS, mean_corr_SIM,mean_corr_pHL, mean_corr_pHN];
% mean_corr_vals_noDC = [mean_corr_pWF,  mean_corr_OSS, mean_corr_pHL];
% mean_corr_vals_DC = [ mean_corr_pWN, mean_corr_SIM, mean_corr_pHN];
% 
% std_pWF_CC = std(pWF_Corrs);
% std_OSS_CC = std(OSS_Corrs);
% std_pHL_CC = std(pHL_Corrs);
% std_pWN_CC = std(pWN_Corrs);
% std_SIM_CC = std(SIM_Corrs);
% std_pHN_CC = std(pHN_Corrs);
% 
% 
% errs_CC = [std_pWF_CC, std_pWN_CC,  std_OSS_CC,  std_SIM_CC, std_pHL_CC, std_pHN_CC]/2;
% errs_CC_noDC = [std_pWF_CC,  std_OSS_CC,   std_pHL_CC]/2;
% errs_CC_DC = [ std_pWN_CC,    std_SIM_CC, std_pHN_CC]/2;
% 
% % b = bar(x,mean_corr_vals, 'facecolor', 'flat');
% % b.CData= clr;
% 
% % b = bar(x2,mean_corr_vals_noDC, 'facecolor', 'flat');
% % b.CData= clr_noDC;
% % 
% b = bar(x2,mean_corr_vals_DC, 'facecolor', 'flat');
% b.CData= clr_noDC;
% 
% 
% 
% yticks([0 0.3 0.6])
% xticks([])
% hold on
% 
% er = errorbar(x2,mean_corr_vals_DC,errs_CC_DC);    
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% 
% 
% hold off
% 
% %%
% 
% 
% 
% % % dFF_pWF = [dFF_pWF1,dFF_pWF2,dFF_pWF3,dFF_pWF4,dFF_pWF5,dFF_pWF6,dFF_pWF7,dFF_pWF8,dFF_pWF9,dFF_pWF10,dFF_pWF11];
% % maxes_OS = zeros(1,11);
% % maxes_pWF = zeros(1,11);
% % % 
% % % for ii = 1:11
% % %     maxes_OS(1,ii) = max(dFF_SIM(:,ii));
% % %     maxes_pWF(1,ii) = max(dFF_pWF(:,ii));
% % % end
% % % avg_max_SIM = mean(maxes_OS);
% % % avg_max_pWF = mean(maxes_pWF);
% % % dFF_ratio = avg_max_SIM/avg_max_pWF;
% % % 
% % % 
% % % 
num_Cells = 6;
fr = 100;
offset_S = 1.5;
WF_scale = 5.0;
HL_scale = 0.55;
r_scale = 0.3;
%xi = 120000; % Paper!
%xf = 360000;
xi = 240000;
xf = 280000;

wi = 47000/10;
wf = 49000/10;
yi = 0.5;
yf = 6.9;
clf
for i = 1:6

    subplot(2,2,4)
    if i <= 10
        plot(xVals,(smooth(dFF_SIM(i,:),3)+i*offset_S)/offset_S,LineWidth= 1.0)
    else
        x3 = linspace(wi,wf,100);
        y3 = ones(100,1)*10;
        plot(x3, y3,'Color','magenta',LineWidth= 5.0)
    end 
    % plot(xVals, dFF_SIM(:,i)+i*offset_S-1,LineWidth= 1.0)
    xlim([xi,xf])
    ylim([yi,yf])
    %xticks([])
    %yticks([])
    %ylabel("SIM")
    %xlabel("Time (ms)")
    ylabel("Cell #")
    hold on
    title("OS-SIM N2N")

    subplot(2,2,1)
    plot(xVals,(smooth(WF_scale*dFF_pWF(i,:),3)+i*offset_S)/offset_S,LineWidth= 1.0)
    %plot(xVals, 2*dFF(dFF_n2n(10,ts:tf)',xVals)+i*offset_S)
    xlim([xi,xf])
    ylim([yi,yf])
    %xticks([])
    %yticks([])
    ylabel("Cell #")
    xlabel("Time (ms)")
    hold on
    title("5.0 * pWF")

    % subplot(2,3,3)
    % plot(xVals,(smooth(HL_scale*dFF_pHL(i,:),3)+i*offset_S)/offset_S,LineWidth= 1.0)
    % %plot(xVals, 2*dFF(dFF_n2n(10,ts:tf)',xVals)+i*offset_S)
    % xlim([xi,xf])
    % ylim([yi,yf])
    % %xticks([])
    % %yticks([])
    % ylabel("Cell #")
    % xlabel("Time (ms)")
    % hold on
    % title("0.55 * pHiLo")

    % subplot(2,3,6)
    % plot(xVals,(smooth(HL_scale*dFF_pHN(i,:),3)+i*offset_S)/offset_S,LineWidth= 1.0)
    % %plot(xVals, 2*dFF(dFF_n2n(10,ts:tf)',xVals)+i*offset_S)
    % xlim([xi,xf])
    % ylim([yi,yf])
    % %xticks([])
    % %yticks([])
    % ylabel("Cell #")
    % xlabel("Time (ms)")
    % hold on
    % title("0.55 * pHiLo N2N")

    subplot(2,2,2)
    plot(xVals,(smooth(2.0*dFF_OSS(i,:),3)+i*offset_S)/offset_S,LineWidth= 1.0)
    %plot(xVals, 2*dFF(dFF_n2n(10,ts:tf)',xVals)+i*offset_S)
    xlim([xi,xf])
    ylim([yi,yf])
    %xticks([])
    %yticks([])
    ylabel("Cell #")
    xlabel("Time (ms)")
    hold on
    title("2.0 * OS-SIM")

    subplot(2,2,3)
    plot(xVals,(smooth(2.5*WF_scale*dFF_pWN(i,:),3)+i*offset_S)/offset_S,LineWidth= 1.0)
    %plot(xVals, 2*dFF(dFF_n2n(10,ts:tf)',xVals)+i*offset_S)
    xlim([xi,xf])
    ylim([yi,yf])
    %xticks([])
    %yticks([])
    ylabel("Cell #")
    xlabel("Time (ms)")
    hold on
    title("7.5 * pWF N2N")


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


% %%
% % % 
% num_Cells =6;
% clf
% corr_SIM = corrcoef(dFF_SIM(1:num_Cells,:)');
% corr_pWF = corrcoef(dFF_pWF(1:num_Cells,:)');
% 
% %corr_WG = corrcoef(dFF_WG(1:num_Cells,:)');
% 
% corr_OSS = corrcoef(dFF_OSS(1:num_Cells,:)');
% corr_pWN = corrcoef(dFF_pWN(1:num_Cells,:)');
% corr_pHL = corrcoef(dFF_pHL(1:num_Cells,:)');
% corr_pHN = corrcoef(dFF_pHN(1:num_Cells,:)');
% 
% 
% 
% mean_corr_SIM = mean(corr_SIM, "all"); 
% mean_corr_pWF = mean(corr_pWF, "all"); 
% mean_corr_OSS = mean(corr_OSS, "all"); 
% mean_corr_pWN = mean(corr_pWN, "all"); 
% mean_corr_pHL = mean(corr_pHL, "all"); 
% mean_corr_pHN = mean(corr_pHN, "all"); 
% 
% 
% 
% %mean_corr_WG = mean(corr_WG, "all");
% %ratio = mean_corr_OS/mean_corr_WF;   
% 
% SIM_Corrs = [];
% pWF_Corrs = [];
% OSS_Corrs = [];
% pWN_Corrs = [];
% pHL_Corrs = [];
% pHN_Corrs = [];
% 
% 
% 
% %hm = zeros(num_Cells,num_Cells);
% %hm2 = zeros(num_Cells,num_Cells);
% 
% for i = 1:num_Cells
%     for l = 1:num_Cells
%         if i < l 
%             %hm(i,l) = corr_OS(i,l);
%             %hm2(i,l) = corr_OS(i,l);
%             SIM_Corrs(end+1) = corr_SIM(i,l);
%             pWF_Corrs(end+1) = corr_pWF(i,l);
%             OSS_Corrs(end+1) = corr_OSS(i,l);
%             pWN_Corrs(end+1) = corr_pWN(i,l);
%             pHL_Corrs(end+1) = corr_pHL(i,l);
%             pHN_Corrs(end+1) = corr_pHN(i,l);
% 
% 
% 
% 
% 
%         elseif i>l
%             %hm(i,l) = corr_WF(i,l);
%             %hm2(i,l) = corr_WG(i,l);
%             %wf_Corrs(end+1) = corr_WF(i,l);
%         else 
%             %hm(i,l) = NaN;
%             %hm2(i,l) = NaN;
%         end
%     end
% end
% 
% % subplot(2,3,1)
% % histogram(pWF_Corrs)
% % title("pWF Corr. Coefs. (Mean = 0.5399)")
% % xlabel("Correlation Coeff.")
% % ylabel("Occurances")
% % 
% % subplot(2,3,2)
% % histogram(OSS_Corrs)
% % title("OS-SIM Corr. Coefs. (Mean = 0.2482)")
% % xlabel("Correlation Coeff.")
% % ylabel("Occurances")
% % 
% % subplot(2,3,3)
% % histogram(pHL_Corrs)
% % title("pHiLo Corr. Coefs. (Mean = 0.1926)")
% % xlabel("Correlation Coeff.")
% % ylabel("Occurances")
% % 
% % subplot(2,3,4)
% % histogram(pWN_Corrs)
% % title("pWF-N2N Corr. Coefs. (Mean = 0.5329)")
% % xlabel("Correlation Coeff.")
% % ylabel("Occurances")
% % 
% % subplot(2,3,5)
% % histogram(SIM_Corrs)
% % title("OS-SIM-N2N Corr. Coefs. (Mean = 0.2357)")
% % xlabel("Correlation Coeff.")
% % ylabel("Occurances")
% % 
% % subplot(2,3,6)
% % histogram(pHN_Corrs)
% % title("pHiLo-N2N Corr. Coefs. (Mean = 0.1887)")
% % xlabel("Correlation Coeff.")
% % ylabel("Occurances")
% % 
% % 
% % 
% % 
% % 
% % 
% 
% 
% 
% 
% 
% 
% % subplot(2,4,3)
% % % % %heatmap(corr_OS,'ColorLimits',[0.1 0.85], 'Colormap', jet)
% % % % %heatmap(corr_OS,'ColorLimits',[0.0 1])
% % heatmap(hm,'ColorLimits',[0.0 1])
% % ylabel("Cell # pWF")
% % xlabel("Cell # pWF")
% % title("Heatmap with pWF")
% % 
% % subplot(2,4,7)
% % % % %heatmap(corr_OS,'ColorLimits',[0.1 0.85], 'Colormap', jet)
% % % % %heatmap(corr_OS,'ColorLimits',[0.0 1])
% % heatmap(hm2,'ColorLimits',[0.0 1])
% % ylabel("Cell # pWF+WGN")
% % xlabel("Cell # pWF+WGN")
% % title("Heatmap with pWF+WGN")
% % 
% % 
% % subplot(2,4,4)
% % y = [mean_corr_OS, mean_corr_WF, mean_corr_WG];
% % x = ["OS" "WF" "WF+WG"];
% % bar(x,y)
% % title("Mean Correlation Coefficient")
% % 
% % subplot(2,4,8)
% % y = [PNR_AVG_OS, PNR_AVG_WF, PNR_AVG_WG ];
% % x = ["OS" "WF" "WF+WG"];
% % bar(x,y)
% % title("Mean Spike Peak-to-Noise Ratio")
% 
% 
% %plot(xVals, dFF_SIM(:,10))
% %xlim([46000,48000])
% % % % % 
% % % % 
% % % % 
% % % % 
% % % % % 
% % subplot(3,3,3)
% % hist(wf_Corrs)
% % xlabel("Corr. Coeff. (N2N)")
% % %xlim([0,0.35])
% % title("Noise2Noise AVG. Coeff = 0.2357")
% % 
% % subplot(3,3,6)
% % hist(os_Corrs)
% % xlabel("Corr. Coeff. (FB)")
% % %xlim([0,0.35])
% % title("Fish Brains AVG. Coeff = 0.2668")
% 
% 
% %% HITS!!
% clf
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
% for z = 12000:35994
%     for i = 1:6
%         if dFF_pWF(i,z) >= max(dFF_pWF(i,12000:35994))*thresh_scale
%             j = j + 1;
%         end
%         hits_pWF(z,1) = j;
% 
%         if dFF_OSS(i,z) >= max(dFF_OSS(i,12000:35994))*thresh_scale
%             k = k + 1;
%         end
%         hits_OSS(z,1) = k;
% 
%         if dFF_pHL(i,z) >= max(dFF_pHL(i,12000:35994))*thresh_scale
%             l = l + 1;
%         end
%         hits_pHL(z,1) = l;
% 
% 
% 
%         if dFF_pWN(i,z) >= max(dFF_pWN(i,12000:35994))*thresh_scale
%             jj = jj + 1;
%         end
%         hits_pWN(z,1) = jj;
% 
%         if dFF_SIM(i,z) >= max(dFF_SIM(i,12000:35994))*thresh_scale
%             kk = kk + 1;
%         end
%         hits_SIM(z,1) = kk;
% 
%         if dFF_pHN(i,z) >= max(dFF_pHN(i,12000:35994))*thresh_scale
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
% 
% 
% frames = linspace(1,35994-12000,35995-12000);
% 
% 
% 
% %subplot(1,3,3)
% % plot(frames, hits_pWF(12000:35994,1),'r',LineWidth=2)
% % hold on 
% % plot(frames,  hits_OSS(12000:35994,1),'b',LineWidth=2)
% % plot(frames, hits_pHL(12000:35994,1),'g',LineWidth=2)
% 
% plot(frames, hits_pWN(12000:35994,1),'r',LineWidth=2)
% hold on 
% plot(frames, hits_SIM(12000:35994,1),'b',LineWidth=2)
% plot(frames, hits_pHN(12000:35994,1),'g',LineWidth=2)
% 
% 
% 
% 
% xlabel("Frames")
% legend("pWF",'OS-SIM','pHiLo', Location = 'northwest')
% %xlim([12000,36000])
% %ylim([yi,yf])
% xticks([10000 20000])
% yticks([3000 6000 9000])
% %ylim([0, max(hits_pWF(:,1))])
% percent = 100*k/j;
% ylabel('Frames above 30 % Threshold')
% % 
% %%
% percent_OS = 100*k/(6*29501);
% percent_WF = 100*j/(6*29501);%%
% percent_HL = 100*l/(6*29501);%%
% %%
% 
% % % t = Tiff('n2n_pHiLo.tif', 'r');
% % % imageData = read_file(t);
% % 
% % filename_hil = 'n2n_pHiLo_16b.tif';
% % [hilData, ~, ~, ~] = read_file(filename_hil);
% % 
% % 
% % %%
% % 
% % 
% % filename_SIM = 'NEIL_1_P1_mc_16s_E_05_Iter_6254_outp_i_interleaved_NOmatchSIM_Reconstruction.tif';
% % filename_pWF = 'NEIL_1_P1_m_interleaved_interleaved_pWF_Reconstruction.tif';
% % filename_hil = 'n2n_pHiLo_16b.tif';
% % 
% % t = Tiff('n2n_pHiLo.tif', 'r');
% % imageData = read(t);
% % 
% % 
% % info = imfinfo(filename_pWF);
% % numframe = 35994;
% % for K = 1 : numframe
% %     rawframes_SIM(:,:,K) = imread(filename_SIM, K);
% %     rawframes_pWF(:,:,K) = imread(filename_pWF, K);
% %     %rawframes_hil(:,:,K) = imread(filename_hil, K);
% % 
% %     %rawframes_rSIM(:,:,K) = imread(filename_rSIM, K);
% % 
% % end
% % simframes = mat2gray(rawframes_SIM);
% % 
% % 
% % pWFframes = mat2gray(rawframes_pWF);
% % %rsimframes = mat2gray(rawframes_rSIM);
% 
% % hilframes = mat2gray(hilData);
% % % 
% % % corr_image_SIM = zeros(247,295);
% % % corr_image_pWF = zeros(247,295);
% % %corr_image_rSIM = zeros(84,230);
% % corr_image_hil = zeros(247,295);
% % 
% % 
% % for x=1:295
% %     for y = 1:247
% %         %cc_sim = corrcoef(simframes(y,x,:),SIMFullROI.Mean);
% %         %cc_pWF = corrcoef(pWFframes(y,x,:),WFFullROI.Mean);
% %         %cc_rsim = corrcoef(rsimframes(y,x,:),rsim(1,:));
% %         cc_hil = corrcoef(hilframes(y,x,:),hilFullROI.Mean);
% % 
% %         %corr_image_sim(y,x) = cc_sim(2,1);
% %         %corr_image_pwf(y,x) = cc_pWF(2,1);
% %         %corr_image_rsim(y,x) = cc_rsim(2,1);
% %         corr_image_hil(y,x) = cc_hil(2,1);
% % 
% % 
% %     end
% % end
% 
% % clims = [0.0 1];
% % subplot(1,3,1)
% % imagesc(corr_image_pwf, clims)
% % colorbar
% % title("pWF Correlation Map")
% % xticks([])
% % yticks([])
% % 
% % subplot(1,3,2)
% % imagesc(corr_image_sim, clims)
% % colorbar
% % xticks([])
% % yticks([])
% % title("OS-SIM Correlation Map")
% % 
% % subplot(1,3,3)
% % imagesc(corr_image_hil, clims)
% % colorbar
% % xticks([])
% % yticks([])
% % title("HiLo Correlation Map")
% 
% % subplot(4,5,15)
% % imagesc(corr_image_rsim)
% % colorbar
% % title("rSIM Correlation Map")
% % 
% % subplot(4,5,20)
% % imagesc(corr_image_hil)
% % colorbar
% % title("pHilo Correlation Map")
% 
% %%  Try Stuff 0402
% clf
% % i = 2
% % p = prctile(smooth(dFF_OSS(i,12000:35994),3),99.9);
% % o = ones(size(dFF_OSS(i,:)));
% % p_arr = p*o;
% % 
% % subplot(2,1,1)
% % plot(xVals(12000:35994),smooth(dFF_OSS(i,12000:35994),3))
% % hold on 
% % plot(xVals(12000:35994),p_arr(12000:35994))
% % 
% % 
% % p = prctile(smooth(dFF_SIM(i,12000:35994),3),99.9);
% % p_arr = p*o;
% % 
% % subplot(2,1,2)
% % plot(xVals(12000:35994),smooth(dFF_SIM(i,12000:35994),3))
% % hold on 
% % plot(xVals(12000:35994),p_arr(12000:35994))
% 
% reconstruction = ["WF" "OS" "HL" "WF-DC" "OS-DC" "HL-DC"];
% rec = table(pnr_pWF',pnr_OSS', pnr_pHL' ,pnr_pWN', pnr_SIM', pnr_pHN', VariableNames = reconstruction);
% matrix_pWF = [pnr_pWF';pnr_pWN'];
% matrix_SIM = [pnr_OSS';pnr_SIM'];
% matrix_pHL = [pnr_pHL';pnr_pHN'];
% matrix_rec = [matrix_pWF, matrix_SIM,matrix_pHL];
% 
% method = [repmat("Widefield",12,1);repmat("OS-SIM",12,1);repmat("HiLo",12,1)];
% denoising = [repmat("None",6,1);repmat("DeepCAD",6,1);repmat("None",6,1);repmat("DeepCAD",6,1);repmat("None",6,1);repmat("DeepCAD",6,1)];
% tbl = table(method,denoising,matrix_rec(:),VariableNames=["ReconstructionMethod" "Denoising" "PeaktoNoiseRatio"]);
% aovLinear = anova(tbl,"PeaktoNoiseRatio ~ ReconstructionMethod + Denoising");
% aovInteraction = anova(tbl,"PeaktoNoiseRatio ~ ReconstructionMethod + Denoising + ReconstructionMethod:Denoising")
% 
% 
% %aov = anova(rec);
% 
% 
% 
% 


%%

% clf
% % y = sim10.Mean(ts:tf,1);
% opol = 15;
% [p,s,mu] = polyfit(xVals',y(:,1),opol);
% f_y = polyval(p,xVals',[],mu);
% 
% dt_ecgnl = y - f_y;
% %dFF = dt_ecgnl./f_y;

perc(8369, 3468)



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
                   'WindowSize', 9000, ...
                   'StepSize', 9000, ...
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




function percent_increase = perc(x,denoised_x)

    percent_increase = (denoised_x-x)/x;

end




