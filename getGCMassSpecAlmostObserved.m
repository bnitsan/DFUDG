function massScaleIn = getGCMassSpecAlmostObserved()

massScaleIn=ReadDataIntoMassSpec(); 

% break most massive GC into smaller ones
heavyGC=max(massScaleIn); % find heavist 
idheavy=find(massScaleIn==heavyGC);
massScaleIn(idheavy)=[];
massScaleIn = [massScaleIn' heavyGC/2 heavyGC/4 heavyGC/4]';

% add three more GCs at smallest mass bin 
% to account for out-of-cut-radius GCs
massScaleIn = [massScaleIn' 1 1.75 2.5]';
end

