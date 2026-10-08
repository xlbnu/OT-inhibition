function get2_2StatisticsValue(xx1,x_num)
xx1_tmp=array2table(xx1,'VariableNames',{'s1','s2','p1','p2'});
% wf={'lambda';'drug'};wl=[2;2];
rm=fitrm(xx1_tmp,'s1,s2,p1,p2 ~ 1',WithinDesign=table([1 -1 1 -1]',...
    [1 1 -1 -1]','VariableNames', {'drug';'context'}));
% mauchly(rm)
ranovatbl = ranova(rm,'WithinModel','context*drug');
effectSizeTable=computeEffectSize(ranovatbl);
hp=pair_ttest(xx1(:,[1 3]),xx1(:,[2 4]));
xx_diff=mean([xx1(:,[1 3])-xx1(:,[2 4]),xx1(:,1)-xx1(:,3)-(xx1(:,2)-xx1(:,4))]);
xx_sd=std([xx1(:,[1 3])-xx1(:,[2 4]),xx1(:,1)-xx1(:,3)-(xx1(:,2)-xx1(:,4))]);
cohen_d=xx_diff./xx_sd;
hp1=pair_ttest(xx1(:,1)-xx1(:,3),(xx1(:,2)-xx1(:,4)));

hp_0=hp;
hp_0(:,2)=hp_0(:,2).*x_num;
hp_0(hp_0(:,2)>1,2)=1;
hp1_0=hp1;
hp1_0(:,2)=hp1_0(:,2).*x_num;
hp1_0(hp1_0(:,2)>1,2)=1;

p_re=[ranovatbl.pValue(3),ranovatbl.pValue(5),ranovatbl.pValue(7)].*x_num;
p_re(p_re>1)=1;

fprintf(['\nsocial : OT-PL pair-test: t(%d)=%.4f, p_adj=%.4f(p=%.4f), Cohen’s d=%.4f; nonsocial : OT-PL pair-test: t(%d)=%.4f, p_adj=%.4f(p=%.4f), Cohen’s d=%.4f;' ...
    ' social-nonsocial : OT-PL pair-test: t(%d)=%.4f, p_adj=%.4f(p=%.4f), Cohen’s d=%.4f;  '],...
    round(hp(1,4)), hp(1,3),hp_0(1,2), hp(1,2),cohen_d(1), round(hp(2,4)), hp(2,3),hp_0(2,2), hp(2,2),cohen_d(2),round(hp1(1,4)), hp1(1,3),hp1_0(1,2), hp1(1,2),cohen_d(3));
fprintf(['drug main effect: F(%d,%d)=%.4f, p_adj=%.4f(p=%.4f), partical η2=%.4f；context main effect: F(%d,%d)=%.4f, p_adj=%.4f(p=%.4f), partical η2=%.4f；' ...
    ' drug×context interaction: F(%d,%d)=%.4f, p_adj=%.4f(p=%.4f), partical η2=%.4f\n\n'],...
    ranovatbl.DF(3),ranovatbl.DF(4),ranovatbl.F(3),p_re(1),ranovatbl.pValue(3),effectSizeTable.PartialEtaSq(2),...
    ranovatbl.DF(5),ranovatbl.DF(6),ranovatbl.F(5),p_re(2),ranovatbl.pValue(5),effectSizeTable.PartialEtaSq(3),...
    ranovatbl.DF(7),ranovatbl.DF(8),ranovatbl.F(7),p_re(3),ranovatbl.pValue(7),effectSizeTable.PartialEtaSq(4));
end


