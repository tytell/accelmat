function plot_bg14_accel_data

flowspeed0 = [1.0 1.0    1.0    1.0    1.0    1.0    1.0    1.5 ...
    1.5    1.5    1.5    1.5    1.5    1.5    1.5    1.5 ...
    1.5    1.75    1.75    1.75    1.75    2.0    2.0    2.0    2.0 ...
    2.0    2.0    2.0    1.75    1.75    1.75    1.75    0.5    0.5];

load('rawdata/bg14/bg14data-analysis.mat');

goodchan = [1 1 1 1 0 1 0 1]' > 0;

good = true;
fn = fieldnames(out);
for i = 1:length(out)
    for j = 1:length(fn)
        if size(out(i).(fn{j}),1) == 6
            switch class(out(i).(fn{j}))
                case 'double'
                    a = NaN(8,size(out(i).(fn{j}),2));
                    
                case 'char'
                    a = repmat('L',8,size(out(i).(fn{j}),2));
                    
                case 'cell'
                    a = cell(8,size(out(i).(fn{j}),2));
            end
            a(goodchan,:) = out(i).(fn{j});
            out(i).(fn{j}) = a;
            good = false;
        elseif isempty(out(i).(fn{j}))
            switch class(out(1).(fn{j}))
                case 'double'
                    out(i).(fn{j}) = NaN;
                case 'char'
                    out(i).(fn{j}) = 'L';
                case 'cell'
                    out(i).(fn{j}) = {[]};
            end
            good = false;
        end
    end
end

if ~good
    save('rawdata/bg14/bg14data-analysis.mat','out');
end

for j = 1:length(fn)
    S.(fn{j}) = catuneven(3, out.(fn{j}));
end
flowspeed = NaN(size(S.speed));
for i = 1:length(out)
    flowspeed(:,:,i) = flowspeed0(i);
end

figureseries('Accel mag');
good = isfinite(S.speed(1,:,:)) & (S.tailbeatnum(1,:,:) <= 5);
plotgroups(flatten(S.tailbeatnum(1,good)), flatten(S.accfwdpk(1,good)), ...
    flatten(flowspeed(1,good)), {'cmf'}, 'means',true, 'legend');
xlabel('Peak accel (g)');

figureseries('Accel mn mag');
good = isfinite(S.speed(1,:,:)) & (S.tailbeatnum(1,:,:) <= 5);
plotgroups(flatten(S.tailbeatnum(1,good)), flatten(S.accfwdmn(1,good)), ...
    flatten(flowspeed(1,good)), {'cmf'}, 'means',true, 'legend');
xlabel('Mean accel (g)');

%------------------------------------------
figureseries('Body amplitude');
colormap(lines(5));
good = isfinite(S.speed(1,:,:));

h(1) = subplot(2,1,1);
plotgroups(S.accfwdpk(1,good), abs(S.tailamp(1,good)), ...
    {flowspeed(1,good)}, {'cm'}, ...
    'means',true, 'legend','regression',true);

h(2) = subplot(2,1,2);
plotgroups(S.accfwdpk(1,good), abs(S.headamp(1,good)), ...
    {flowspeed(1,good)}, {'cm'}, ...
    'means',true, 'legend','regression',true);

linkaxes(h,'xy');
axis([-0.02 0.4 0 0.18]);

%------------------------------------------
figureseries('Body amplitude - accmn');
colormap(lines(5));
good = isfinite(S.speed(1,:,:));

h(1) = subplot(2,1,1);
plotgroups(S.accfwdmn(1,good), abs(S.tailamp(1,good)), ...
    {flowspeed(1,good)}, {'cm'}, ...
    'means',true, 'legend','regression',true);

h(2) = subplot(2,1,2);
plotgroups(S.accfwdmn(1,good), abs(S.headamp(1,good)), ...
    {flowspeed(1,good)}, {'cm'}, ...
    'means',true, 'legend','regression',true);

linkaxes(h,'xy');
axis([-0.02 0.4 0 0.18]);

%------------------------------------------
figureseries('Tail freq');
colormap(lines(5));
good = isfinite(S.speed(1,:,:));

plotgroups(S.accfwdpk(1,good), abs(S.freq(1,good)), ...
    {flowspeed(1,good)}, {'cm'}, ...
    'means',true, 'legend','regression',true);
axis([-0.02 0.4 2 9]);

%------------------------------------------
figureseries('Tail freq - accmn');
colormap(lines(5));
good = isfinite(S.speed(1,:,:));

plotgroups(S.accfwdmn(1,good), abs(S.freq(1,good)), ...
    {flowspeed(1,good)}, {'cm'}, ...
    'means',true, 'legend','regression',true);
axis([-0.02 0.4 2 9]);

%------------------------------------------
figureseries('Burst amp');
colormap(lines(4));
showchan = find(goodchan);
for i = 1:6
    hburstamp(i) = subplot(3,2,i);
    
    c = showchan(i);
    good = isfinite(S.burstamp(c,:)) & ...
        (flowspeed(c,:) > 0.5) & (S.accfwdpk(c,:) < 0.4);
    plotgroups(S.accfwdpk(c,good), S.burstamp(c,good), ...
        {flowspeed(c,good)},{'cmf'}, 'legend','regression',true);
end
linkaxes(hburstamp,'xy');

%------------------------------------------
figureseries('Burst amp3');
colormap(lines(4));
c = 3;
good = isfinite(S.burstamp(c,:)) & ...
    (flowspeed(c,:) > 0.5) & (S.accfwdpk(c,:) < 0.4);
plotgroups(S.accfwdpk(c,good), S.burstamp(c,good), ...
    {flowspeed(c,good)},{'cmf'}, 'legend','regression',true);


%------------------------------------------
figureseries('Burst amp3 - accmn');
colormap(lines(4));
c = 3;
good = isfinite(S.burstamp(c,:)) & ...
    (flowspeed(c,:) > 0.5) & (S.accfwdpk(c,:) < 0.4);
plotgroups(S.accfwdmn(c,good), S.burstamp(c,good), ...
    {flowspeed(c,good)},{'cmf'}, 'legend','regression',true);


%------------------------------------------
figureseries('Burst duty');
colormap(lines(4));
showchan = find(goodchan);
for i = 1:6
    hburstduty(i) = subplot(3,2,i);
    
    c = showchan(i);
    good = isfinite(S.burstduty(c,:)) & ...
        (flowspeed(c,:) > 0.5) & (S.accfwdpk(c,:) < 0.4);
    plotgroups(S.accfwdpk(c,good), S.burstduty(c,good), ...
        {flowspeed(c,good)},{'cmf'}, 'legend','regression',true);
end
linkaxes(hburstduty,'xy');

%------------------------------------------
figureseries('Burst duty 1 6');
colormap(lines(4));
showchan = [1 6];
for i = 1:2
    hburstduty16(i) = subplot(2,1,i);
    
    c = showchan(i);
    good = isfinite(S.burstduty(c,:)) & ...
        (flowspeed(c,:) > 0.5) & (S.accfwdpk(c,:) < 0.4);
    plotgroups(S.accfwdpk(c,good), S.burstduty(c,good), ...
        {flowspeed(c,good)},{'cmf'}, 'legend','regression',true);
end
linkaxes(hburstduty16,'xy');

%------------------------------------------
figureseries('Burst duty 1 6 - accmn');
colormap(lines(4));
showchan = [1 6];
for i = 1:2
    hburstduty16(i) = subplot(2,1,i);
    
    c = showchan(i);
    good = isfinite(S.burstduty(c,:)) & ...
        (flowspeed(c,:) > 0.5) & (S.accfwdpk(c,:) < 0.4);
    plotgroups(S.accfwdmn(c,good), S.burstduty(c,good), ...
        {flowspeed(c,good)},{'cmf'}, 'legend','regression',true);
end
linkaxes(hburstduty16,'xy');

%------------------------------------------
figureseries('Burst on');
colormap(lines(4));
showchan = find(goodchan);
for i = 1:6
    hburston(i) = subplot(3,2,i);
    
    c = showchan(i);
    
    on1 = S.burstonphase(c,:);
    on1 = mod(on1 + 0.6, 1) - 0.6;
    
    good = isfinite(on1) & ...
        (flowspeed(c,:) > 0.5) & (S.accfwdpk(c,:) < 0.4);
    plotgroups(S.accfwdpk(c,good), on1(good), ...
        {flowspeed(c,good)},{'cmf'}, 'legend','regression',true);
end
linkaxes(hburston,'xy');

%------------------------------------------
figureseries('Burst on - accmn');
colormap(lines(4));
showchan = find(goodchan);
for i = 1:6
    hburston(i) = subplot(3,2,i);
    
    c = showchan(i);
    
    on1 = S.burstonphase(c,:);
    on1 = mod(on1 + 0.6, 1) - 0.6;
    
    good = isfinite(on1) & ...
        (flowspeed(c,:) > 0.5) & (S.accfwdpk(c,:) < 0.4);
    plotgroups(S.accfwdmn(c,good), on1(good), ...
        {flowspeed(c,good)},{'cmf'}, 'legend','regression',true);
end
linkaxes(hburston,'xy');

%------------------------------------------
figureseries('report');

good = isfinite(S.tailamp(1,:,:)) & (flowspeed(1,:,:) == 1.75) & ...
    (S.accfwdpk(1,:,:) >= 0);

hrep(1) = subplot(2,2,1);
mplot(S.accfwdpk(1,good), abs(S.tailamp(1,good)), 'ko', 'MarkerSize',5);
fprintf('Amplitude\n');
plotpolybounds(S.accfwdpk(1,good), abs(S.tailamp(1,good)), 'showstats',true);
axis tight;
xlim([0 0.45]);
ylabel('Tail amplitude (L)');

hrep(2) = subplot(2,2,2);
mplot(S.accfwdpk(1,good), S.freq(1,good), 'ko', 'MarkerSize',5);
fprintf('Frequency\n');
plotpolybounds(S.accfwdpk(1,good), S.freq(1,good), 'showstats',true);
axis tight;
xlim([0 0.45]);
ylabel('Tail beat frequency (Hz)');

c = 1;

hrep(3) = subplot(2,2,3);
good = isfinite(S.burstduty(c,:,:)) & (flowspeed(c,:,:) == 1.75) & ...
    (S.accfwdpk(c,:,:) >= 0);

mplot(S.accfwdpk(1,good), S.burstduty(c,good), 'bsf', 'MarkerSize',5);
fprintf('Midbody duty\n');
plotpolybounds(S.accfwdpk(1,good), S.burstduty(1,good), 'color','b', 'showstats',true);
axis tight;
xlim([0 0.45]);
ylabel('Burst duty cycle');

hrep(4) = subplot(2,2,4);

on1 = S.burstonphase(c,:);
on1 = mod(on1 + 0.6, 1) - 0.6;

mplot(S.accfwdpk(1,good), on1(good), 'bsf', 'MarkerSize',5);
fprintf('Tail duty\n');
plotpolybounds(S.accfwdpk(1,good), on1(good), 'color','b', 'showstats',true);
axis tight;
xlim([0 0.45]);
ylabel('Burst onset phase');

c = 6;

subplot(2,2,3);
good = isfinite(S.burstduty(c,:,:)) & (flowspeed(c,:,:) == 1.75) & ...
    (S.accfwdpk(c,:,:) >= 0);

addmplot(S.accfwdpk(c,good), S.burstduty(c,good), 'rx', 'MarkerSize',5);
fprintf('Midbody onset\n');
plotpolybounds(S.accfwdpk(c,good), S.burstduty(c,good), 'color','r', 'showstats',true);
axis tight;
xlim([0 0.45]);
ylabel('Burst duty cycle');
xlabel('Peak acceleration (g)');

subplot(2,2,4);

on1 = S.burstonphase(c,:);
on1 = mod(on1 + 0.6, 1) - 0.6;

addmplot(S.accfwdpk(c,good), on1(good), 'rx', 'MarkerSize',5);
fprintf('Tail onset\n');
plotpolybounds(S.accfwdpk(1,good), on1(good), 'color','r', 'showstats',true);
axis tight;
xlim([0 0.45]);
ylabel('Burst onset phase');
xlabel('Peak acceleration (g)');

set(hrep,'box','off','tickdir','out');






