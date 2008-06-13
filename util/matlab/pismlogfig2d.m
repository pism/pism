function pismlogfig2d(x,y,A,lev);

if nargin > 3
    A = min(A,max(lev));
    A = max(A,min(lev));
end
surf(x,y,flipud(A'),log10(flipud(A')),'LineStyle','none')
axis tight
view(2)
h = colorbar;
if nargin > 3
    levlabels = cell(1,length(lev));
    for k=1:length(lev)
        levlabels{k} = num2str(lev(k));
    end
    set(h,'YTick',log10(lev))
    set(h,'YTickLabel',levlabels)
end