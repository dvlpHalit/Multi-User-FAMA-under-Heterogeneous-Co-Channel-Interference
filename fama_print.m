function fama_print(h,name)
%FAMA_PRINT  Save a figure in both formats the paper build may need.
%   name.pdf  picked up automatically by pdflatex
%   name.eps  for a latex+dvips route, and for journal submission systems
print(h,'-dpdf' ,'-r600',[name '.pdf']);
print(h,'-depsc2','-r600',[name '.eps']);
fprintf('  wrote %s.pdf and %s.eps\n',name,name);
end
