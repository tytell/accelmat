function sync_video_and_data(vidfiles, datafiles, varargin)

opt.videooffset = 0;
opt.trigger = 'end';
opt.triggerfrac = 1.0;
opt.fps = [];
opt.plotvars = {};
opt.isoverlay = true;
opt = parsevarargin(opt, varargin, 3);

if ((nargin == 0) || isempty(vidfiles))
    [fn,pn] = uigetfile('*.avi','Choose video(s)', 'MultiSelect','on');
    if (iscell(fn))
        vidfiles = cell(size(fn));
        for i = 1:length(fn)
            vidfiles{i} = fullfile(pn,fn{i});
        end
    else
        vidfiles = {fullfile(pn,fn)};
    end
end

if ((nargin < 2) || isempty(datafiles))
    [fn,pn] = uigetfile('*.mat','Choose data file(s)', 'MultiSelect','on');
    if (iscell(fn))
        datafiles = cell(size(fn));
        for i = 1:length(fn)
            datafiles{i} = fullfile(pn,fn{i});
        end
    else
        datafiles = {fullfile(pn,fn)};
    end
end

if ischar(vidfiles)
    vidfiles = {vidfiles};
end

nvid = length(vidfiles);
vidnames = cell(1,nvid);

for i = 1:nvid
    vid(i) = VideoReader2(vidfiles{i});
    fps(i) = vid(i).FrameRate;
    nframes(i) = vid(i).NumberOfFrames;

    [~,fn] = fileparts(vidfiles{i});
    vidnames{i} = fn;

    if ~isempty(opt.fps)
        if numel(opt.fps) == 1
            newfps = opt.fps;
        else
            assert(length(opt.fps) == nvid);
            newfps = opt.fps(i);
        end
        
        if ischar(opt.trigger)
            trigtype = opt.trigger;
        elseif iscell(opt.trigger) && (numel(opt.trigger) == 1)
            trigtype = opt.trigger{1};
        else
            assert(length(opt.trigger) == nvid);
            trigtype = opt.trigger{i};
        end
        
        if numel(opt.triggerfrac) == 1
            trigfrac = opt.triggerfrac;
        else
            assert(length(opt.triggerfrac) == nvid);
            trigfrac = opt.triggerfrac(i);
        end            
    else
        prompt = {'Frames per second:', 'Trigger time (fraction):'};
        def = {num2str(fps(i)), '1.0'};
        
        vals = inputdlg(prompt, vidnames{i}, 1, def);
        if isempty(vals)
            fprintf('Cancelling...\n');
            return;
        end
        
        newfps = str2double(vals{1});
        trig = vals{2};
        trig = str2double(vals{2});
        if (isnan(trig) && ischar(vals{2}))
            trigtype = vals{2};
        else
            trigtype = 'frac';
            trigfrac = trig;
        end
    end
    
    fps(i) = newfps;
    
    tvid{i} = (0:nframes(i)-1)/fps(i);
    switch trigtype
        case 'end'
            tvid{i} = tvid{i} - tvid{i}(end);
        case 'start'
             % do nothing
        case {'frac','fraction'}
            assert((trigfrac >= 0) && (trigfrac <= 1));
            ind = round((nframes(i)-1) * trigfrac)+1;
            tvid{i} = tvid{i} - tvid{i}(ind);
        otherwise
            error('Unknown trigger type %s', opt.trigger);
    end
end

if ischar(datafiles)
    datafiles = {datafiles};
end

ndata = length(datafiles);
if ~isempty(opt.plotvars) && iscell(opt.plotvars)
    if ischar(opt.plotvars{1})
        opt.plotvars = {opt.plotvars};
    end
    assert(length(opt.plotvars) == ndata);
end
    
nsubplot = 0;
for i = 1:ndata
    vars = who('-file',datafiles{i});
    [pn,fn] = fileparts(datafiles{i});
    
    if isempty(opt.plotvars)
        ist = strncmpi('t',vars,1);
        vars1 = [vars(ist); vars(~ist)];
        
        [tsel,status] = listdlg('PromptString','Select time variable:', ...
            'SelectionMode','single', ...
            'ListString', vars1, 'InitialValue',1);
        if ~status
            fprintf('Cancelling...\n');
            return;
        end

        tvar = vars1{tsel};

        vars1 = [vars(~ist); vars(ist)];
        [varsel,status] = listdlg('PromptString','Select variables to plot:', ...
            'SelectionMode','multiple', ...
            'ListString', vars1);
        if ~status
            fprintf('Cancelling...\n');
            return;
        end

        plotvars = vars1(varsel);

        if (length(varsel) > 1)
            over = questdlg('Overlay plots?','Plot type','Yes','No','Yes');
            isover = strcmp(over,'Yes');
        else
            isover = true;
        end
    elseif ~isempty(opt.plotvars)
        assert(length(opt.plotvars{i}) >= 2);
        tvar = opt.plotvars{i}{1};
        plotvars = opt.plotvars{i}(2:end);
        isover = opt.isoverlay;
    end

    if isover
        plotcmd1 = cell(1,2*length(plotvars)+1);
        plotcmd1{1} = @plot;
        for j = 1:length(plotvars)
            plotcmd1{2*j} = tvar;
            plotcmd1{2*j + 1} = plotvars{j};
        end
        plotcmds{i}{1} = plotcmd1;
        plotnames{i}{1} = plotvars;
        nsubplot = nsubplot + 1;
    else
        plotcmd1 = cell(1,length(plotvars));
        for j = 1:length(plotvars)
            plotcmd1{j} = {@plot,tvar,plotvars{j}};
        end
        plotcmds{i} = plotcmd1;
        plotnames{i} = plotvars;
        nsubplot = nsubplot + numel(plotcmd1);
    end
    F{i} = load(datafiles{i},tvar,plotvars{:});
end

fig = openfig(mfilename, 'new');
data = guihandles(fig);

data.fig = fig;
data.tvid = tvid;
data.vid = vid;
data.fps = fps;
data.vidfiles = vidfiles;
data.datafiles = datafiles;
data.vidnames = vidnames;
data.htimer = -1;

data = do_layout(data, nvid,nsubplot);

tshow = 0;
trngvid = [-Inf Inf];
for i = 1:nvid
    [~,frload] = min(abs(tvid{i} - tshow));

    im = read(vid(i), frload);
    him(i) = image(im, 'Parent',data.hvidax(i));
    axis(data.hvidax(i),'equal','tight','off','ij');

    lab = sprintf('%s: Frame %d/%d', data.vidnames{i}, frload,length(tvid{i}));
    set(data.hvidlab(i), 'String',lab);
    
    trngvid(1) = max(trngvid(1), min(tvid{i}));
    trngvid(2) = min(trngvid(2), max(tvid{i}));
end
step = 1./max(fps);

plotind = 1;
trngplot = [-Inf Inf];
for i = 1:length(plotcmds)
    for j = 1:length(plotcmds{i})
        plotcmd1 = plotcmds{i}{j};
        for k = 2:length(plotcmd1)
            if isfield(F{i},plotcmd1{k})
                plotcmd1{k} = F{i}.(plotcmd1{k});
            end
        end
        hln{i,j} = feval(plotcmd1{1},data.hplotax(plotind),plotcmd1{2:end});
        set(hln{i,j},'HitTest','off');
        set(data.hplotax(plotind), 'ButtonDownFcn',@on_click_axes);
        ylabel(data.hplotax(plotind), plotnames{i}{j});
        
        axis(data.hplotax(plotind), 'tight');
        hbar(plotind) = vertplot(data.hplotax(plotind), tshow, 'k--');
        
        xl = xlim(data.hplotax(plotind));
        trngplot(1) = max(trngplot(1), xl(1));
        trngplot(2) = min(trngplot(2), xl(2));
        plotind = plotind + 1;
    end
end
linkaxes(data.hplotax, 'x');

trng = [max(trngvid(1),trngplot(1)) min(trngvid(2),trngplot(2))];
set(data.hplotax, 'XLim',trng);

set(data.slider, 'Min',trng(1), 'Max',trng(2), 'Value',tshow, ...
    'SliderStep',[step/(trng(2)-trng(1)) 0.1], ...
    'Callback',@on_slider);
set(data.playbutton, 'Callback', @on_play_clicked);
set(data.startbutton, 'Callback', @on_start_clicked);
set(data.endbutton, 'Callback', @on_end_clicked);
set(data.quitbutton, 'Callback', @on_quit_clicked);

data.him = him;
data.hln = hln;
data.hbar = hbar;
data.trng = trng;
data.tshow = tshow;
data.step = step;

guidata(fig, data);

set(fig,'KeyPressFcn',@on_key_press);

uiwait(fig);

if (ishandle(fig))
    set(fig,'KeyPressFcn','');
    delete(fig);
end


function data = do_layout(data, nvidax,nplotax)

gap = 0.02;

h1 = data.vidaxes1;
hlab1 = data.vidlabel1;
pos0 = get(h1,'Position');
w = (pos0(3) - gap*(nvidax-1)) / nvidax;

hvidax = zeros(1,nvidax);
hvidlab = zeros(1,nvidax);
for i = 1:nvidax
    x = pos0(1) + (i-1)*(w+gap);

    if (i > 1)
        hvidax(i) = copyobj(h1, data.fig);
        hvidlab(i) = copyobj(hlab1, data.fig);
    else
        hvidax(i) = h1;
        hvidlab(i) = hlab1;
    end
    
    pos = pos0;
    pos([1 3]) = [x w];
    set(hvidax(i),'Position',pos);

    
    pos = get(hvidlab(i), 'Position');
    pos([1 3]) = [x w];
    set(hvidlab(i), 'Position',pos, ...
        'String',sprintf('%s (frame %d/%d)', data.vidnames{i}, 1,length(data.tvid{i})));
end
data.hvidax = hvidax;
data.hvidlab = hvidlab;

h1 = data.plotaxes1;
pos0 = get(h1,'Position');
h = (pos0(4) - gap*(nplotax-1)) / nplotax;

hplotax = zeros(1,nplotax);
for i = 1:nplotax
    y = pos0(2)+pos0(4) - (i-1)*(h+gap) - h;
    if (i > 1)
        hplotax(i) = copyobj(h1, data.fig);
    else
        hplotax(i) = h1;
    end
    
    pos = pos0;
    pos([2 4]) = [y h];
    set(hplotax(i),'Position',pos);
end
data.hplotax = hplotax;



function on_key_press(src, event)

data = guidata(src);
switch event.Key
    case 'leftarrow'
        if (data.tshow - data.step >= data.trng(1))
            data.tshow = data.tshow-data.step;
            data = update(data);
        end
        
    case 'rightarrow'
        if (data.tshow + data.step <= data.trng(2))
            data.tshow = data.tshow+data.step;
            data = update(data);
        end
        
    case 'q'
        uiresume(data.fig);
end
guidata(src,data);


function on_click_axes(src,event)

data = guidata(src);
c = get(src, 'CurrentPoint');
data.tshow = c(1,1);

data = update(data);

guidata(src,data);


function on_slider(src,event)

data = guidata(src);

v = get(src,'Value');
data.tshow = v;

data = update(data);

guidata(src,data);


function on_start_clicked(src,event)

data = guidata(src);
data.tshow = data.trng(1);
data = update(data);
guidata(src,data);


function on_end_clicked(src,event)

data = guidata(src);
data.tshow = data.trng(end);
data = update(data);
guidata(src,data);


function on_quit_clicked(src,event)

data = guidata(src);
uiresume(data.fig);


function on_play_clicked(src,event)

data = guidata(src);

rate = get(data.playspeededit, 'String');
rate = str2double(rate);
if (get(src,'Value') > 0)
    if ((rate > 0) && (rate < 1000))
        playstart = clock;
        data.htimer = timer('TimerFcn',{@on_play_timer,data,rate,data.tshow,playstart}, ...
            'Period',1/rate, 'ExecutionMode','fixedRate', ...
            'StopFcn', {@on_stop_timer,data,rate,data.tshow,playstart});
        start(data.htimer);
    else
        set(src,'Value',0);
    end
else
    if (isvalid(data.htimer))
        stop(data.htimer);
        delete(data.htimer);
        data.htimer = -1;
        
        data = guidata(src);
    end
end
guidata(src,data);



function on_play_timer(src,event, data,rate,vidtime0,playtime0)

data.tshow = vidtime0 + etime(event.Data.time, playtime0)*rate/max(data.fps);
if (data.tshow > data.trng(2))
    stop(src);
end
update(data);


function on_stop_timer(src,event, data,rate,vidtime0,playtime0)

data.tshow = vidtime0 + etime(event.Data.time, playtime0)*rate/max(data.fps);
if (data.tshow > data.trng(2))
    data.tshow = data.trng(2);
end
guidata(data.fig, data);


function data = update(data)

data.tshow = diground(data.tshow,0.001);
for i = 1:length(data.vid)
    [~,frload] = min(abs(data.tvid{i} - data.tshow));

    im = read(data.vid(i), frload);
    set(data.him(i), 'CData',im);
    
    set(data.hvidlab(i), 'String',sprintf('%s (frame %d/%d)', data.vidnames{i}, frload,length(data.tvid{i})));    
    set(data.timeedit, 'String',sprintf('%.3f',data.tshow));
end

set(data.hbar, 'XData', [data.tshow; data.tshow]);
set(data.slider, 'Value', data.tshow);

