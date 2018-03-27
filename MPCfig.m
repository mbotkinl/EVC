function MPCfig(data,i)

%H=gca();
%figure
% fig=get(groot,'CurrentFigure');
% if isempty(fig)
%     t1=1;
%     t2=length(data);
% else
%     t1
% end


%change old lines to dotted and faded
fig=get(groot,'CurrentFigure');
if ~isempty(fig)
    H=gca();
    lines=H.Children;
    for ii=1:length(lines)
        lines(ii).LineStyle=':';
        lines(ii).Color(4) = 1-(ii*.1); % older predictions get lighter and lighter (assume 1 is most recent)???
    end
end
plot(1+i:length(data)+i,data,'b')
hold on
