function mdl = wm_train_svm_ecoc(X,y)
    classes = unique(y(:))';
    W = wm_class_weights(y);
    t = templateSVM('KernelFunction','linear', 'BoxConstraint',1, 'Standardize',false);
    mdl = fitcecoc(X, y, 'Learners', t, 'Coding', 'onevsone', ...
        'ClassNames', classes, 'Weights', W);
end
