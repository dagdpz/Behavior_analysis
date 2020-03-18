function [Out Info]=Simple_memory_analysis(MA_in)

Info.sessions=unique([MA_in.selected.session]);

SU=[MA_in.binary.success];

LS=real([MA_in.saccades.tar_pos])<=0;
RS=real([MA_in.saccades.tar_pos])>0;
CH=[MA_in.binary.choice];
IN=~CH;

conditions=[LS&IN&SU; RS&IN&SU; LS&CH&SU; RS&CH&SU];
labels={'LI','RI','LC','RC'};
for l=1:numel(labels)
    idx=conditions(l,:);
    fn=labels{l};
    Out.RT.(fn)=nanmean([MA_in.saccades(idx).lat]);
    Out.N.(fn)=sum(idx);
    Out.ACX.(fn)=nanmean(real([MA_in.saccades(idx).endpos]-[MA_in.saccades(idx).tar_pos]));
    Out.ACY.(fn)=nanmean(imag([MA_in.saccades(idx).endpos]-[MA_in.saccades(idx).tar_pos]));
end
sum_C=Out.N.LC+Out.N.RC;
sum_I=Out.N.LI+Out.N.RI;
 Out.N.LC=Out.N.LC/(sum_C);
 Out.N.RC=Out.N.RC/(sum_C);
 Out.N.LI=Out.N.LI/(sum_I);
 Out.N.RI=Out.N.RI/(sum_I);

FN=fieldnames(Out);
summary_figure1=figure('units','normalized','outerposition',[0 0 1 1]);
for s=1:numel(FN)
subplot(2,2,s)
hold on
for l=1:numel(labels)
    bar(l,Out.(FN{s}).(labels{l}));    
end
title([num2str(Info(1).sessions) '  ' FN{s}])
set(gca,'xtick',1:numel(labels),'xticklabel',labels);

end

end