function sync_video_and_data(vidfile, datafile, varargin)

opt.videooffset = 0;
opt.trigger = 'end';
opt = parsevarargin(opt, varargin, 3);

vid = VideoReader2(vidfile);
fps = vid.FrameRate;
nframes = vid.NumberOfFrames;

newfps = input(sprintf('Frames per second? (default = %d) ',fps));
if ~isempty(newfps)
    fps = newfps;
end

frames = 0:nframes-1;
switch opt.trigger
    case 'end'
        frames = frames - frames(end);
    case 'start'
         % do nothing
    otherwise
        error('Unknown trigger type %s', opt.trigger);
end

vars = who('-file',datafile);
[pn,fn] = fileparts(datafile);
fprintf('Found the following variables in %s:\n', fn);

fprintf('   %s\n', vars{:});

fprintf('\n\n');
if ~ismember('t',vars) || ~inputyn('Use t as time variable?')
    tvar = input('Time variable? ','s');
    
    if ~ismember(tvar,vars)
        error('%s is not a variable in the file.',tname);
    end
else
    tvar = 't';
end

plotvar = input('Plot variable? ','s');

F = load(datafile,tvar,plotvar);

fig = figure;
hax(1) = subplot(3,1,1:2);

frshow = 0;
[~,frload] = min(abs(frames - frshow));

im = read(vid, frload);
him = image(im);
axis equal tight off ij

hax(2) = subplot(3,1,3);
hln = plot(F.(tvar), F.(plotvar));
hbar = vertplot(frshow/fps, 'k--');

data.fig = fig;
data.hax = hax;
data.him = him;
data.hln = hln;
data.hbar = hbar;
data.frames = frames;
data.frshow = frshow;
data.vid = vid;
data.fps = fps;

guidata(fig, data);

set(fig,'KeyPressFcn',@on_key_press);

uiwait(fig);

set(fig,'KeyPressFcn','');


function on_key_press(src, event)

data = guidata(src);
switch event.Key
    case 'leftarrow'
        if (data.frshow > data.frames(1))
            data.frshow = data.frshow-1;
            data = update(data);
        end
        
    case 'rightarrow'
        if (data.frshow < data.frames(end))
            data.frshow = data.frshow+1;
            data = update(data);
        end
        
    case 'q'
        uiresume(data.fig);
end
guidata(src,data);



function data = update(data)

[~,frload] = min(abs(data.frames - data.frshow));

im = read(data.vid, frload);
set(data.him, 'CData',im);

t1 = data.frshow/data.fps;

set(data.hbar, 'XData', [t1; t1]);
