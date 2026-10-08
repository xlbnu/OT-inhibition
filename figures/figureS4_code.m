%%%%%%%%%%%%%%%%%%% figure S4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;
% Locate bundled inputs relative to this saved script, independent of pwd.
figureDir = fileparts(mfilename('fullpath'));
codeRoot = fileparts(figureDir);
dataRoot = fullfile(codeRoot,'figure_data');
addpath(fullfile(codeRoot,'utilities'));
f_file = [fullfile(codeRoot,'figure_outputs','figureS4') filesep];
if ~isfolder(f_file)
    mkdir(f_file);
end
load(fullfile(dataRoot,'model_results','twolambda_model_fit_parameter_85.mat'),'ots_minb1','pls_minb1','otn_minb1','pln_minb1');

ots_xx1=cat(1,ots_minb1.optParams);
pls_xx1=cat(1,pls_minb1.optParams);
otn_xx1=cat(1,otn_minb1.optParams);
pln_xx1=cat(1,pln_minb1.optParams);
%
pn=[7 8 1:6 9 10];p_params=nan(10,1);
for i=1:10
[~,p_params(i,1)]=ttest(ots_xx1(:,pn(i))-otn_xx1(:,pn(i))-(pls_xx1(:,pn(i))-pln_xx1(:,pn(i))));
end

%% figure S4B
% two lambda
p_ylim=[-0 0 0 0 0 0 0 0 0;1.5 1.5 1.5 15 25 0.2 0.2 1500 1500]';
paramNames0 = {'baseline drift', 'baseline self bias', 'baseline other bias', 'self input gain', 'partner input gain', 'noise strength', 'mutual inhibition', 'mutual inhibition', 'urgency time constant', 'urgency onset'};
paramNames = {'v_0', 'b_s/b_0', 'b_p/b_0', 'k_s','k_p', 'σ', 'λ_o_→_s', 'λ_s_→_o', '\tau', 't_0'};
paramNames2 = {'v_0', 'b_s', 'b_o', 'k1', 'k2','σ', 'λ_(o→s)', 'λ_s→o', 'τ', 't_0'};
pn=[7 8 1:6 9 10];
for i=1:10
    ax1=figureX1([23 20]);
    hold on;
    x_pos=figureTmp1_OT([ots_xx1(:,pn(i)),pls_xx1(:,pn(i)),otn_xx1(:,pn(i)),pln_xx1(:,pn(i))],10); 
    setFixedPlotArea(ax1, [15 15]);
    hold off;box off;set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'}); gcaf1(gca,35);
    ylabel([paramNames0{pn(i)},' (',paramNames{pn(i)},')'],'FontSize',40,'FontName','Arial');
     setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
    print(gcf,[f_file,'TwoLambda-modelParameter(',paramNames2{(pn(i))},').svg'], '-dsvg','-vector','-r600');
end


