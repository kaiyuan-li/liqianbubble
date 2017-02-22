function figformat()
% figure format

                set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 24);
                set(gca,'XMinorTick','on')
                set(gca,'YMinorTick','on')
                %set(gca,'YTick',[0:1:5])
                set(gcf,'Color',[1 1 1])
                set(gcf, 'units', 'inches', 'position', [2 2 8 4])
                %set(gca,'units,'inches','position', [2 2 10 15])
a=gca;
a = a.Children;
for ii = 1:length(a)
    set(a(ii),'linewidth',2);
end            
end

