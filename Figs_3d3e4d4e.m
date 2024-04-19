%% Figure 3d, 3e, 4d, 4e

% Step 1, drag the .csv files from '\\contrast' into workspace

% Step 2, 
clf 

contrast = @(x) (x/(max(x)));
calc_contrast = @(x) (x(13)-x(1))/(x(13)+x(1));

[W, BW] = adjust_baseline(pWFnoDC.Gray_Value);
[WN, BWN] = adjust_baseline(pWF.Gray_Value);
[O, BO] = adjust_baseline(OSnoDC.Gray_Value);
ON = adjust_baseline(OS.Gray_Value);
H = adjust_baseline(pHiLonoDC.Gray_Value);
HN = adjust_baseline(pHiLo.Gray_Value);

%% Fig 3D
clf
plot(contrast(W)-min(contrast(W)),'r', LineWidth=2);
hold on 
plot(contrast(O)-min(contrast(O)),'b', LineWidth=2);
plot(contrast(H)-min(contrast(H)),'g', LineWidth=2);
yticks([0  0.5 1.5])
xticks([0 100 200])
%xlim([1 26])
legend(["WF" "OS" "HL"], NumColumns=2)

%% Fig 4D
clf 
plot(contrast(WN)-min(contrast(WN)),'r', LineWidth=2)
hold on 
plot(contrast(ON)-min(contrast(ON)),'b', LineWidth=2)
plot(contrast(HN)-min(contrast(HN)),'g', LineWidth=2)

%% Fig 3E
clf
CC_WF = calc_contrast(W);
CC_OS = calc_contrast(O);
CC_HL = calc_contrast(H);

CC_noDC = [CC_WF  CC_OS  CC_HL ];
clr_noDC = [1 0 0; 0 0 1;  0 1 0 ];

b = bar(CC_noDC, "FaceColor","flat");
b.CData= clr_noDC;

xticks([])
yticks([0.5])

%% Fig 4E
clf
CC_WN = calc_contrast(WN);
CC_ON = calc_contrast(ON);
CC_HN = calc_contrast(HN);

CC_DC = [CC_WN  CC_ON  CC_HN ];

clr_noDC = [1 0 0; 0 0 1;  0 1 0 ];

b = bar(CC_DC, "FaceColor","flat");
b.CData= clr_noDC;

xticks([])
yticks([0.5])


%% Adjust baseline function

function [corrected, Baseline] = adjust_baseline(x)
    corrected = linspace(1,209,209);
    WF_og2 = x(1:209);
    WF11 = WF_og2(1,1);
    WF2091 = WF_og2(end);
    for i = 1:length(WF_og2)
        Baseline(i) = WF11 + i*(WF11-WF2091)/209;
    end
    for i = 1:209
        corrected(i) = WF_og2(i)-Baseline(210-i)/2;
    end 

end