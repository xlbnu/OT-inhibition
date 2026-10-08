function report = check_figure_inputs()
% CHECK_FIGURE_INPUTS Read-only checks of the bundled figure input files.
%   Run after setup_figure_paths. No figures or analyses are executed.
%   Verifies paths and explicitly requested top-level MAT variables, plus
%   BrainNet surfaces/configurations and S7 sliding-window MAT/NIfTI assets.
%   Does not validate scientific results, nested fields, or panel rendering.
codeRoot = fileparts(mfilename('fullpath'));
requests = jsondecode(fileread(fullfile(codeRoot,'figures','figure_inputs.json')));
report = struct('script',{},'file',{},'ok',{},'message',{});
for k = 1:numel(requests)
    r = requests(k);
    p = fullfile(codeRoot,'figure_data',r.folder,r.file);
    entry = struct('script',r.script,'file',p,'ok',false,'message','');
    try
        assert(isfile(p),'Input file is missing.');
        info = whos('-file',p);
        missing = setdiff(string(r.variables),string({info.name}));
        assert(isempty(missing),'Missing variables: %s',strjoin(missing,', '));
        entry.ok = true;
        entry.message = 'File and requested variables are available.';
    catch ME
        entry.message = ME.message;
    end
    report(end+1) = entry; %#ok<AGROW>
end
assets = jsondecode(fileread(fullfile(codeRoot,'figures','figure_assets.json')));
for k = 1:numel(assets)
    r = assets(k);
    p = fullfile(codeRoot,strrep(r.path,'/',filesep));
    entry = struct('script',r.script,'file',p,'ok',false,'message','');
    try
        assert(isfile(p),'Required figure asset is missing.');
        if ~isempty(r.variables)
            info = whos('-file',p);
            missing = setdiff(string(r.variables),string({info.name}));
            assert(isempty(missing),'Missing variables: %s',strjoin(missing,', '));
        end
        entry.ok = true;
        entry.message = 'Figure asset and requested variables are available.';
    catch ME
        entry.message = ME.message;
    end
    report(end+1) = entry; %#ok<AGROW>
end
fprintf('Figure resource checks passed: %d/%d\n',sum([report.ok]),numel(report));
if any(~[report.ok])
    disp(struct2table(report(~[report.ok])));
    warning('OT_code:FigureInputs','Resolve the listed input problems before plotting.');
end
end
