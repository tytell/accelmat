function plotPrelimEMGacc

vidfile = 'rawdata/2014-05-14/Ventral video/Bg3_17_0Hz_V_012.avi';
kinfile = 'rawdata/2014-05-14/Ventral video/Bg3_17_0Hz_V_012.mat';
accfile = 'rawdata/2014-05-14/bg3_17_0Hz_012acc.mat';
emgfile = 'rawdata/2014-05-14/bg3_17_0Hz_012emg.mat';

finbeatframe = [2549 2624 2673 2740 2786 2829 2899 2962 3010 3091];
showbeats = [-11 -10];

imax = [920 1620 286 641];

gmmss = 9800;       % mm/s^2
tshow = [-3 0];

kin = load(kinfile,'tymm','humms','t','haxmmss');
kin.t0 = kin.t(end);
kin.t = kin.t - kin.t0;

acc = load(accfile,'t','acchi','orient');

emg = load(emgfile,'t','LTail','RTail');

[~,indleft] = findpeaks2(kin.tymm,'max');
tailleftt = kin.t(indleft);

[~,indright] = findpeaks2(kin.tymm,'min');
tailrightt = kin.t(indright);

showframes = [indleft(end+showbeats); indright(end+showbeats)];
showframes = showframes(:);

vid = VideoReader2(vidfile);
I = zeros(vid.Height,vid.Width, length(showframes));
for i = 1:length(showframes);
    I1 = read(vid,showframes(i));
    I(:,:,i) = im2double(I1(:,:,1));
end

clf;
n = length(showframes);
for i = 1:n
    him(i) = subplot(8,n, i+[0 n 2*n]);
    imshow(I(imax(3):imax(4),imax(1):imax(2),i),'InitialMagnification','fit');
    addplot(100,100,'r*');
end

tailleftt = tailleftt((tailleftt >= tshow(1)) & (tailleftt <= tshow(2)));
tailrightt = tailrightt((tailrightt >= tshow(1)) & (tailrightt <= tshow(2)));

iskin = (kin.t >= tshow(1)) & (kin.t <= tshow(end));
isacc = (acc.t >= tshow(1)) & (acc.t <= tshow(end));
isemg = (emg.t >= tshow(1)) & (emg.t <= tshow(2));

h(1) = subplot(8,1,4);
plot(kin.t(iskin),150-kin.tymm(iskin),'k');
axis tight;
ylabel({'Tail position','(mm)'});
xtick labeloff;

h(2) = overlayplot('top');
addplot(h(2),acc.t(isacc),acc.orient(isacc,3)*180/pi,'b--');
axis tight;
vertplot(h(2),tailleftt,'k--', tailrightt,'k:');
yl = ylim;
addplot(kin.t(showframes),yl(2),'r*');
ylabel({'Yaw angle','(deg)'});
xtick labeloff;

h(3) = subplot(8,1,5);
plot(acc.t(isacc),acc.acchi(isacc,2),'r-');
addplot(kin.t(iskin),-kin.haxmmss(iskin) / gmmss);
axis tight;
vertplot(tailleftt,'k--', tailrightt,'k:');
ylabel({'Forward','acceleration (g)'});
xtick labeloff;

h(4) = subplot(8,1,6);
plot(emg.t(isemg), emg.LTail(isemg));
axis tight;
vertplot(tailleftt,'k--', tailrightt,'k:');
ytick off;
xtick labeloff;

h(5) = subplot(8,1,7);
plot(emg.t(isemg),emg.RTail(isemg));
axis tight;
vertplot(tailleftt,'k--', tailrightt,'k:');
ytick off;
xtick labeloff;

h(6) = subplot(8,1,8);
plot(acc.t(isacc),acc.acchi(isacc,1),'g-');
axis tight;
vertplot(kin.t(finbeatframe),'g-');
addscalebars('x','xlen',0.2,'units','sec');
ylabel({'Vertical','acceleration (g)'});

linkaxes(h,'x');
xlim(tshow);
set(h,'Box','off');

