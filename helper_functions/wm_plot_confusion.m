function wm_plot_confusion(cm, labels, outPng, titleStr)
    fig = figure('Color','w','Position',[100 100 650 580]);
    imagesc(cm);
    axis square;
    colorbar;
    set(gca,'XTick',1:numel(labels),'XTickLabel',labels);
    set(gca,'YTick',1:numel(labels),'YTickLabel',labels);
    xlabel('Predicted label');
    ylabel('True label');
    title(titleStr, 'Interpreter','none');

    for r = 1:size(cm,1)
        for c = 1:size(cm,2)
            val = cm(r,c);
            if val > 0.55
                col = 'w';
            else
                col = 'k';
            end
            text(c, r, sprintf('%.2f', val), ...
                'HorizontalAlignment','center', 'FontWeight','bold', 'Color',col);
        end
    end

    saveas(fig, outPng);
    savefig(fig, strrep(outPng,'.png','.fig'));
    close(fig);
end
