% COntrast PLots 031924

clf

contrast = @(x) (x/(max(x)));
calc_contrast = @(x) (x(13)-x(1))/(x(13)+x(1));



%%  0402 analysis
clf

[pWF_adj, pWF_bl] = adjust_baseline2(pWFnoDC.Gray_Value);
[OSS_adj, OSS_bl] = adjust_baseline2(OSnoDC.Gray_Value);
[HiLo_adj, HiLo_bl] = adjust_baseline2(pHiLonoDC.Gray_Value);
[pWF_adj_dc, pWF_bl_dc] = adjust_baseline2(pWF.Gray_Value);
[OSS_adj_dc, OSS_bl_dc] = adjust_baseline2(OS.Gray_Value);
[HiLo_adj_dc, HiLo_bl_dc] = adjust_baseline2(pHiLo.Gray_Value);

% plot(pWF_adj, 'r')
% hold on 
% plot(OSS_adj, 'b')
% plot(HiLo_adj, 'g')


plot(pWF_adj_dc, 'r')
hold on 
plot(OSS_adj_dc, 'b')
plot(HiLo_adj_dc, 'g')


% 
% 
% plot(pWF.Gray_Value, 'r')
% hold on
% plot(OS.Gray_Value, 'b')
% plot(pHiLo.Gray_Value, 'g')


xlim([0 27])


legend(["pWF" "OS" "pHiLo"])



%%
% 
% [W, BW] = adjust_baseline(pWFnoDC.Gray_Value);
% [WN, BWN] = adjust_baseline(pWF.Gray_Value);
% [O, BO] = adjust_baseline(OSnoDC.Gray_Value);
% ON = adjust_baseline(OS.Gray_Value);
% H = adjust_baseline(pHiLonoDC.Gray_Value);
% HN = adjust_baseline(pHiLo.Gray_Value);






% 

% WF_Correct = linspace(1,209,209);
% WF_og = WFnoDC.Gray_Value;
% WF1 = WF_og(1,1);
% WF209 = WF_og(209);
% for i = 1:209
%     WF_Baseline(i) = WF1 + i*(WF1-WF209)/209;
%     WF_Corrected(i) = (WF_Baseline(i) + WF_og(i))/2;
% end 
% plot(WF_Baseline)
% hold on 
% plot(WFnoDC.Gray_Value,'r', LineWidth=2)
% plot(WF_Corrected)

% 
% % 
% % %subplot(1,2,1)
% plot(contrast(W)-min(contrast(W)),'r', LineWidth=2);
% hold on 
% %plot(contrast(WN)-min(contrast(WN)),'r', LineWidth=2)
% %hold on
% plot(contrast(O)-min(contrast(O)),'b', LineWidth=2);
% %plot(contrast(ON)-min(contrast(ON)),'b', LineWidth=2)
% plot(contrast(H)-min(contrast(H)),'g', LineWidth=2);
% %plot(contrast(HN)-min(contrast(HN)),'g', LineWidth=2)
% % 
% % plot(256*(WF.Gray_Value-(WF.Gray_Value(1))),'r--', LineWidth=2);
% % hold on 
% % plot(256*(WFnoDC.Gray_Value-(WFnoDC.Gray_Value(1))),'r', LineWidth=2)
% % plot(256*(OS.Gray_Value-(OS.Gray_Value(1))),'b--', LineWidth=2);
% % plot(OSnoDC.Gray_Value-(OSnoDC.Gray_Value(1)),'b', LineWidth=2)
% % plot(pHiLo.Gray_Value-(pHiLo.Gray_Value(1)),'g--', LineWidth=2);
% % plot(pHiLonoDC.Gray_Value-(pHiLonoDC.Gray_Value(1)),'g', LineWidth=2)
% % 
% % plot(W,'r--', LineWidth=2);
% % hold on 
% % plot(WN,'r', LineWidth=2)
% % plot(O,'b--', LineWidth=2);
% % plot(ON,'b', LineWidth=2)
% % % plot(H,'g--', LineWidth=2);
% % % plot(HN,'g', LineWidth=2)
% 
% yticks([0  0.5 1.5])
% xticks([0 100 200])
% xlim([1 26])
% legend(["WF" "OS" "HL"], NumColumns=2)

% % 

% subplot(1,4,3)
% CC_WF = calc_contrast(W);
% CC_WN = calc_contrast(WN);
% CC_OS = calc_contrast(O);
% CC_ON = calc_contrast(ON);
% CC_HL = calc_contrast(H);
% CC_HN = calc_contrast(HN);
% 
% CC = [CC_WF CC_WN CC_OS CC_ON CC_HL CC_HN];
% CC_noDC = [CC_WF  CC_OS  CC_HL ];
% CC_DC = [CC_WN  CC_ON  CC_HN ];
% 
% x2 = [1 2 3];
% clr_noDC = [1 0 0; 0 0 1;  0 1 0 ];
% 
% b = bar(CC_DC, "FaceColor","flat");
% b.CData= clr_noDC;
% 
% xticks([])
% yticks([0.5])

function [corrected, Baseline] = adjust_baseline(x)
    corrected = linspace(1,209,209);
    WF_og2 = x(1:209);
    WF11 = WF_og2(1,1);
    WF2091 = WF_og2(end);
    for i = 1:length(WF_og2)
        Baseline(i) = WF11 + i*(WF11-WF2091)/209;
    end
    for i = 1:209
        %corrected = (smooth(x,50) - min(smooth(x,50)))/(max(smooth(x,50))+min(smooth(x,50)));
        corrected(i) = WF_og2(i)-Baseline(210-i)/2;
    end 

end

function [corrected, Baseline] = adjust_baseline2(x)
    corrected = linspace(1,26,26);
    og = x(1:26);
    og1 = og(1);
    og_end = og(end);
    for i = 1:length(og)
        Baseline(i) = og1 - i*(og1-og_end)/26;
    end
    for i = 1:26
        corrected(i) = (og(i)-Baseline(i))/(max(og)+min(og));
    end 

end




% 

%%%%%%%%%%%%%%%%%%%%%%%%%%% 200Hz 
% COntrast PLots 031924

% clf
% 
% contrast = @(x) (x/(max(x)));
% calc_contrast = @(x) (x(13)-x(1))/(x(1)+x(13));
% 
% 
% [W, BW] = adjust_baseline(WFnoDC.Gray_Value);
% [WN, BWN] = adjust_baseline(WF.Gray_Value);
% [O, BO] = adjust_baseline(OSnoDC.Gray_Value);
% ON = adjust_baseline(OS.Gray_Value);
% H = adjust_baseline(pHiLonoDC.Gray_Value);
% HN = adjust_baseline(pHiLo.Gray_Value);
% 
% 
% 
% 
% subplot(1,2,1)
% plot(contrast(W)-min(contrast(W)),'r--', LineWidth=2);
% hold on 
% plot(contrast(WN)-min(contrast(WN)),'r', LineWidth=2)
% plot(contrast(O)-min(contrast(O)),'b--', LineWidth=2);
% plot(contrast(ON)-min(contrast(ON)),'b', LineWidth=2)
% plot(contrast(H)-min(contrast(H)),'g--', LineWidth=2);
% plot(contrast(HN)-min(contrast(HN)),'g', LineWidth=2)
% 
% 
% yticks([0  0.5 1.5])
% xticks([0 100 200])
% xlim([1 209])
% legend(["WF" "WN" "OS" "ON" "HL" "HN"], NumColumns=2)
% 
% % 
% 
% subplot(1,4,3)
% CC_WF = calc_contrast(W);
% CC_WN = calc_contrast(WN);
% CC_OS = calc_contrast(O);
% CC_ON = calc_contrast(ON);
% CC_HL = calc_contrast(H);
% CC_HN = calc_contrast(HN);
% 
% CC = [CC_WF CC_WN CC_OS CC_ON CC_HL CC_HN];
% 
% bar(CC)
% xticks([])
% yticks([0.5])
% 
% function [corrected, Baseline] = adjust_baseline(x)
%     corrected = linspace(1,209,209);
%     WF_og2 = x(1:209);
%     WF11 = WF_og2(1,1);
%     WF2091 = WF_og2(end);
%     for i = 1:length(WF_og2)
%         Baseline(i) = WF11 + i*(WF11-WF2091)/209;
%     end
%     for i = 1:209
%         %corrected = (smooth(x,50) - min(smooth(x,50)))/(max(smooth(x,50))+min(smooth(x,50)));
%         corrected(i) = WF_og2(i)-Baseline(210-i)/2;
%     end 
% 
% end
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%% 500Hz 
% COntrast PLots 031924

% clf
% 
% contrast = @(x) (x/(max(x)));
% calc_contrast = @(x) (x(127)-x(137))/(x(137)+x(127));
% 
% 
% [W, BW] = adjust_baseline(WFnoDC.Gray_Value,137);
% [WN, BWN] = adjust_baseline(WF.Gray_Value,137);
% [O, BO] = adjust_baseline(OSnoDC.Gray_Value,137);
% ON = adjust_baseline(OSDC.Gray_Value,137);
% H = adjust_baseline(pHiLonoDC.Gray_Value,137);
% HN = adjust_baseline(pHiLo.Gray_Value,137);
% 
% 
% subplot(1,2,1)
% plot(contrast(W)-min(contrast(W)),'r--', LineWidth=2);
% hold on 
% plot(contrast(WN)-min(contrast(WN)),'r', LineWidth=2)
% plot(contrast(O)-min(contrast(O)),'b--', LineWidth=2);
% plot(contrast(ON)-min(contrast(ON)),'b', LineWidth=2)
% plot(contrast(H)-min(contrast(H)),'g--', LineWidth=2);
% plot(contrast(HN)-min(contrast(HN)),'g', LineWidth=2)
% 
% yticks([0  0.5 1.5])
% xticks([0 50 100])
% xlim([1 137])
% legend(["WF" "WN" "OS" "ON" "HL" "HN"], NumColumns=2)
% 
% % 
% 
% subplot(1,4,3)
% CC_WF = calc_contrast(W);
% CC_WN = calc_contrast(WN);
% CC_OS = calc_contrast(O);
% CC_ON = calc_contrast(ON);
% CC_HL = calc_contrast(H);
% CC_HN = calc_contrast(HN);
% 
% CC = [CC_WF CC_WN CC_OS CC_ON CC_HL CC_HN];
% 
% bar(CC)
% xticks([])
% yticks([0.5])
% 
% function [corrected, Baseline] = adjust_baseline(x,p)
%     corrected = linspace(1,p,p);
%     WF_og2 = x(1:p);
%     WF11 = WF_og2(1,1);
%     WF2091 = WF_og2(end);
%     for i = 1:length(WF_og2)
%         Baseline(i) = WF11 + i*(WF11-WF2091)/p;
%     end
%     for i = 1:p
%         %corrected = (smooth(x,50) - min(smooth(x,50)))/(max(smooth(x,50))+min(smooth(x,50)));
%         corrected(i) = WF_og2(i)-Baseline(p+1-i)/2;
%     end 
% 
% end
% 
% 
% 

