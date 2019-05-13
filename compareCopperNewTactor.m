%% compare the new digital touch probe with the copper tape short circuit
%
% 8.10.2018 David.J.Caldwell

%load('G:\My Drive\buttonVSTact/FSR_characterize_tapeDown.mat')

%files = {'G:\My Drive\copperTapeRTvalidate\8_10_2018_newTactor\responseTiming-2.mat',...,
%};

% tactor vs FSR taped down

[file,path] = uigetfile('-file');
%%
count = 1;
fsrPeaksLocsTotal = [];
load(fullfile(path,file))

coppernData = Tact.data(:,1);
tactData = Tact.data(:,2);

%%
fsTact = Tact.info.SamplingRateHz;

t = 1e3*(0:length(tactData)-1)/fsTact;

figure
plot(t,tactData,'linewidth',2)
hold on
plot(t,coppernData,'linewidth',2)
legend({'digital touch probe','electrical short circuit'})
ylabel('Signal (V)')
xlabel('time (ms)')
set(gca,'fontsize',14);
title('Electrical short circuit vs. digital touch probe')

%%

tButton = (0:size(tactData,1)-1)/fsTact;
tactData(tactData >= 0.015) = 2; % threshold at 0.015
[tactPks,tactLocs] = findpeaks(tactData,tButton,'MinpeakDistance',0.4,'Minpeakheight',1.9);
[tactPksSamps,tactLocsSamps] = findpeaks(tactData,'MinpeakDistance',fsTact*0.4,'Minpeakheight',1.9);

figure
findpeaks(tactData,tButton,'MinpeakDistance',0.4,'Minpeakheight',1.9);

%%

% additional parameters
postStim = 100;
sampsPostStim = round(postStim/1e3*fsTact);

preStim = 100;
sampsPreStim = round(preStim/1e3*fsTact);
tEpoch = 1e3*round([-sampsPreStim:sampsPostStim-1])/fsTact;

tactorEpoched = squeeze(getEpochSignal(tactData,tactLocsSamps-sampsPreStim,tactLocsSamps+ sampsPostStim));
fsrEpoched = squeeze(getEpochSignal(coppernData,tactLocsSamps-sampsPreStim,tactLocsSamps+ sampsPostStim));
figure
subplot(2,1,1)
plot(tEpoch,tactorEpoched)
subplot(2,1,2)
plot(tEpoch,fsrEpoched)
%%
% take diff of fsrEpoched

diffFsrEpoched = [diff(fsrEpoched); zeros(1,size(fsrEpoched,2)) ];

diffFsrEpoched(diffFsrEpoched>0) = 0;
diffFsrEpoched = abs(diffFsrEpoched);

diffFsrEpoched(diffFsrEpoched>=0.15) = 0.15;
diffFsrEpoched(diffFsrEpoched<0.15) = 0;

% loop through the fsr epoched, find the transition point as it starts to
% go negative

fsrPeakLocs = [];
for trial = 1:size(diffFsrEpoched,2)
    [pks,locs] = findpeaks(diffFsrEpoched(:,trial),tEpoch);
    
    if ~isempty(locs)
        fsrPeakLocs(trial) = locs(1);
    else
        fsrPeakLocs(trial) = nan;
    end
    clearvars pks locs
end
%%
figure
histogram(fsrPeakLocs,'numBins',15)
title('FSR vs Tactor')
ylabel('count')
xlabel('difference between the analog FSR and digital tactor (ms)')
set(gca,'fontsize',14);
nanmedian(fsrPeakLocs)

fsrPeaksLocsTotal = [fsrPeaksLocsTotal fsrPeakLocs];
count = count + 1;
%%

fsrPeaksLocsTotal = fsrPeaksLocsTotal(fsrPeaksLocsTotal>-40 & fsrPeaksLocsTotal<0);
figure
histogram(fsrPeaksLocsTotal,'numBins',30)
title('Electrical short circuit vs. digital touch probe')
ylabel('count')
xlabel({'difference between the electrical short circuit',' and digital touch probe onset (ms)'})
set(gca,'fontsize',14);
nanmedian(fsrPeaksLocsTotal)
vline(nanmedian(fsrPeaksLocsTotal),'r',['median = ' num2str(nanmedian(fsrPeaksLocsTotal)) ' ms'])
vline(nanmean(fsrPeaksLocsTotal),'m',['mean = ' num2str(nanmean(fsrPeaksLocsTotal)) ' ms'])

nanmedian(fsrPeaksLocsTotal)
nanmean(fsrPeaksLocsTotal)
std(fsrPeaksLocsTotal((~isnan(fsrPeaksLocsTotal)) & (fsrPeaksLocsTotal<0)))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return

