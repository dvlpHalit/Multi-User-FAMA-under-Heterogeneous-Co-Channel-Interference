function fama_style(h,xl,yl,lg,loc)
%FAMA_STYLE  Uniform IEEE-ready styling. ASCII-only labels, LaTeX interpreter.
set(h,'Units','inches','Position',[1 1 3.45 2.75]);   % single column
set(gca,'FontSize',8,'FontName','Times','Box','on','LineWidth',0.6);
xlabel(xl,'Interpreter','latex','FontSize',9);
ylabel(yl,'Interpreter','latex','FontSize',9);
if nargin>3 && ~isempty(lg)
    if nargin<5, loc='best'; end
    legend(lg,'Interpreter','latex','FontSize',7,'Location',loc,'Box','on');
end
grid on; set(gca,'GridAlpha',0.25,'MinorGridAlpha',0.12);
end
