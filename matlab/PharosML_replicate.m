clear all % 

PharosML_initial
PharosML_setup_for_NN
PharosML_train_NN
NNresult_check



return

function [] = RecordTime() 
persistent TimeRec
    if isempty(TimeRec)
        TimeRec = now;
    else
        TimeRec = [TimeRec,now];
    end
end