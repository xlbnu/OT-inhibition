function result = png_frames_to_gif(frameFiles,outputFile,frameRate)
% PNG_FRAMES_TO_GIF Assemble exported S7 RGB PNG frames chronologically.
%   result = png_frames_to_gif(renderResult.frame_files,outputFile,1);
%   A frame directory may also be supplied instead of an explicit file list.
%   Filenames must encode window times: frame_01_t_m0p80_m0p60_*.png.
%   Existing PNGs are not modified. A common palette is used for all frames.
if nargin<3, frameRate = 1; end
validateattributes(frameRate,{'numeric'},{'scalar','finite','positive','<=',100});
if ischar(frameFiles) || (isstring(frameFiles) && isscalar(frameFiles))
    assert(isfolder(frameFiles),'Expected a frame directory or file list.');
    files = dir(fullfile(frameFiles,'frame_*.png'));
    frameFiles = fullfile({files.folder},{files.name});
elseif isstring(frameFiles)
    frameFiles = cellstr(frameFiles);
end
assert(iscellstr(frameFiles) && ~isempty(frameFiles),'No PNG frames supplied.');
[outputDir,~,extension] = fileparts(outputFile);
assert(strcmpi(extension,'.gif'),'Output filename must end in .gif.');
windows = zeros(numel(frameFiles),2);
for k = 1:numel(frameFiles)
    assert(isfile(frameFiles{k}),'Missing PNG: %s',frameFiles{k});
    [~,base] = fileparts(frameFiles{k});
    tokens = regexp(base,'^frame_\d+_t_([mp]?\d+p\d+)_([mp]?\d+p\d+)(?:_|$)', ...
        'tokens','once');
    assert(~isempty(tokens),'Unrecognized frame time: %s',base);
    for j = 1:2
        token = tokens{j};
        if startsWith(token,'p'), token = token(2:end); end
        windows(k,j) = str2double(strrep(strrep(token,'m','-'),'p','.'));
    end
end
[~,order] = sort(windows(:,1));
windows = windows(order,:); frameFiles = frameFiles(order);
assert(all(diff(windows(:,1))>0) && all(windows(:,2)>windows(:,1)), ...
    'Frame windows must have unique starts and positive durations.');
frames = cell(numel(frameFiles),1);
for k = 1:numel(frameFiles)
    [rgb,map,alpha] = imread(frameFiles{k});
    assert(isempty(map) && isempty(alpha) && isa(rgb,'uint8') && size(rgb,3)==3, ...
        'Expected an RGB PNG exported by the S7 rendering functions.');
    if k>1
        assert(isequal(size(rgb),size(frames{1})),'All PNG frame dimensions must match.');
    end
    frames{k} = rgb;
end
samples = cellfun(@(rgb) rgb(round(linspace(1,size(rgb,1),128)), ...
    round(linspace(1,size(rgb,2),128)),:),frames,'UniformOutput',false);
[~,palette] = rgb2ind(cat(1,samples{:}),256,'nodither');
if ~isempty(outputDir) && ~isfolder(outputDir), mkdir(outputDir); end
for k = 1:numel(frames)
    indexed = rgb2ind(frames{k},palette,'nodither');
    if k==1
        imwrite(indexed,palette,outputFile,'gif','LoopCount',Inf,'DelayTime',1/frameRate);
    else
        imwrite(indexed,palette,outputFile,'gif','WriteMode','append','DelayTime',1/frameRate);
    end
end
result = struct('output_file',outputFile,'frame_count',numel(frames), ...
    'frame_files',{frameFiles},'time_windows',windows,'frame_rate',frameRate);
fprintf('GIF saved: %s (%d chronologically ordered frames)\n',outputFile,numel(frames));
end
