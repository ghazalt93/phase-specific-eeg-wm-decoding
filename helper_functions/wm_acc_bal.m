function [acc,bal,rec,cm] = wm_acc_bal(y,yh,classes)
    cm = confusionmat(y,yh,'Order',classes);
    acc = sum(diag(cm))/max(sum(cm(:)),1);
    rec = diag(cm)./max(sum(cm,2),1);
    bal = mean(rec,'omitnan');
end
