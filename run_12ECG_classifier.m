function [score, label] = run_12ECG_classifier(data,header_data,classes, model)

    num_classes = length(classes);

    label = zeros([1,num_classes]);
    score = ones([1,num_classes]);
    
    
    
       % Use your classifier here to obtain a label and score for each class.
% % %     features = get_12ECG_features(data,header_data);
% % %     
% % %     score = mnrval(model,features);		
% % %     [~,idx] = max (score);
% % %     
% % %     label(idx)=1;
    
   % fprintf(' before get_12ECG_cls_ivo size(data)=%6.0f%6.0f\n');
    
    [score, label]=get_12ECG_cls_ivo(data);
    
 
    
    
    
    
end



