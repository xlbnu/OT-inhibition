%%%%%%%%%%%%%%%%%%% figure S5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;
% Locate bundled inputs relative to this saved script, independent of pwd.
figureDir = fileparts(mfilename('fullpath'));
codeRoot = fileparts(figureDir);
dataRoot = fullfile(codeRoot,'figure_data');
addpath(fullfile(codeRoot,'utilities'));
f_file = [fullfile(codeRoot,'figure_outputs','figureS5') filesep];
if ~isfolder(f_file)
    mkdir(f_file);
end
load(fullfile(dataRoot,'MEG_results','trajectory_sourceEFR_GLM_beta_ROIgamma.mat'),'sb_chosen0','sb_chosen','cb_chosen0','cb_chosen',...
    'otsbs_chosen','otsbc_chosen','mb_chosen','otsb_chosen1','plsb_chosen1','otnb_chosen1','plnb_chosen1','vchan_time');


%% figure S5A
% brainNet to generate gamma-band source ∩ oxytocin reciptor gene
% expression map
cfg = struct();
cfg.brainnet_dir = fullfile(codeRoot,'utilities','BrainNet-Viewer');
cfg.surface_file = fullfile(cfg.brainnet_dir,'BrainMesh_ICBM152Left.nv');
cfg.cfg_file = fullfile(dataRoot,'MEG_results','BrainNet_cfg_medial_left.mat');

cfg.mapping_dir = fullfile(dataRoot,'MEG_results');
cfg.output_dir = f_file;

% Render the map assigned to this panel.
cfg.input_pattern = 'decision_gamma_band_source_OXTRconjunction_.nii';

% Preserve saved display thresholds and color limits.
cfg.fix_global_color_scale = false;
cfg.save_tif_only = true;
result = batch_brainnet_nii_to_tif(cfg);

%% figure S5B
int_color1 = [187 181 214]./255;
int_color  = [132 121 181]./255;

grayLight = [224 218 208]./255;
grayDark  = [ 78  69  59]./255;

stay_sem_color   = [0.88 0.85 0.94];
change_sem_color = [0.88 0.85 0.94];
mixed_sem_color  = [0.86 0.86 0.86];

% =========================================================
% Single-participant signed beta: subject × time
% ==========================================================
stay_self_sub = squeeze(otsbs_chosen(:,:,2));
stay_other_sub = squeeze(otsbs_chosen(:,:,3));

change_self_sub = squeeze(otsbc_chosen(:,:,2));
change_other_sub = squeeze(otsbc_chosen(:,:,3));

mixed_self_sub = squeeze(mb_chosen(:,:,2));
mixed_other_sub = squeeze(mb_chosen(:,:,3));

% time window
stay_idx = vchan_time>-0.7 & vchan_time<-0.35;
change_idx = vchan_time>-0.4 & vchan_time<-0.10;
mixed_idx = vchan_time>-0.7 & vchan_time<-0.0;

% =========================================================
% Bootstrap trajectory
% ==========================================================
stay_result = bootstrap_trajectory_sem(stay_self_sub,stay_other_sub,stay_idx, ...
    'NBoot',5000, 'FirstSmoothMethod','movmean', 'FirstSmoothWindow',20, ...
    'SecondSmoothMethod','gaussian','SecondSmoothWindow',20,'Seed',2026);

change_result = bootstrap_trajectory_sem(change_self_sub,change_other_sub,change_idx, ...
    'NBoot',5000,'FirstSmoothMethod','movmean','FirstSmoothWindow',20, ...
    'SecondSmoothMethod','movmean','SecondSmoothWindow',20,'Seed',2027);

mixed_result = bootstrap_trajectory_sem(mixed_self_sub,mixed_other_sub,mixed_idx, ...
    'NBoot',5000,'FirstSmoothMethod','movmean','FirstSmoothWindow',20, ...
    'SecondSmoothMethod','movmean','SecondSmoothWindow',20,'Seed',2028);

% =========================================================
% Generate gradient colors correctly
% ==========================================================
fraction_s = linspace(0,1,nnz(stay_idx))';
gradientColors_s = int_color1 + fraction_s.*(int_color-int_color1);

fraction_c = linspace(0,1,nnz(change_idx))';
gradientColors_c = int_color1 + fraction_c.*(int_color-int_color1);

fraction_m = linspace(0,1,nnz(mixed_idx))';
gradientColors_m = grayLight + fraction_m.*(grayDark-grayLight);

% =========================================================
% Plotting
% ==========================================================
ax1 = figureX1([24 21]);
hold(ax1,'on');

% First draw the SEM region, then the mean line
draw_trajectory_sem_union(ax1, mixed_result.upper_x,mixed_result.upper_y, ...
    mixed_result.lower_x,mixed_result.lower_y,mixed_sem_color,0.8,true);

draw_trajectory_sem_union(ax1,stay_result.upper_x,stay_result.upper_y, ...
    stay_result.lower_x,stay_result.lower_y,stay_sem_color,0.35,true);

draw_trajectory_sem_union(ax1,change_result.upper_x,change_result.upper_y, ...
    change_result.lower_x,change_result.lower_y,change_sem_color,0.35,true);

% Gradient-colored mean trajectory
plot_gradient_path(ax1,mixed_result.mean_x,mixed_result.mean_y,gradientColors_m,false,8);

plot_gradient_path(ax1,stay_result.mean_x,stay_result.mean_y,gradientColors_s,false,8);

plot_gradient_path(ax1,change_result.mean_x,change_result.mean_y,gradientColors_c,false,8);

xlabel(ax1,'magnitude(self)');
ylabel(ax1,'magnitude(partner)');

xline(ax1,0,'k--');
yline(ax1,0,'k--');

ylim(ax1,[0 0.07]);
xlim(ax1,[0 0.07]);

% self = other
plot(ax1,[0 0.06],[0 0.06], ...
    'k--','LineWidth',1.5, ...
    'HandleVisibility','off');

yticks(ax1,[0:0.035:0.07]);
xticks(ax1,[0:0.035:0.07]);

axis(ax1,'square');
box(ax1,'off');

setFixedPlotArea(ax1,[15 15]);gcaf1(ax1,35);
setAxesForPPT(ax1,36,2);
hold(ax1,'off');
print(gcf,[f_file,'MEG-source-conjunction-trajectory_sem.svg'], '-dsvg','-vector','-r600');
%% figure S5C
% interaction time course
ot_color1 = [226 106 83;242 199 178; 88 118 227;171 198 251]./255;int_color1=[132 121 181]./255;
diff_stay_social = smoothdata(squeeze((otsb_chosen1(:,:,2))),2,'movmean',20)-smoothdata(squeeze((plsb_chosen1(:,:,2))),2,'movmean',20);
diff_stay_asocial = smoothdata(squeeze((otnb_chosen1(:,:,2))),2,'movmean',20)-smoothdata(squeeze((plnb_chosen1(:,:,2))),2,'movmean',20);
ax1=figureX1([25 21]);
hold on
errorbound(vchan_time,mean(diff_stay_social,1,'omitnan'),std(diff_stay_social,'omitnan')./sqrt(size(diff_stay_social,1)),'color',ot_color1(1,:),'linewidth',5);
errorbound(vchan_time,mean(diff_stay_asocial,1,'omitnan'),std(diff_stay_asocial,'omitnan')./sqrt(size(diff_stay_asocial,1)),'color',ot_color1(3,:),'linewidth',5);
plot_column_significance(vchan_time,diff_stay_social-diff_stay_asocial,'Color',int_color1,'YOffsetRatio',0.4,'LineWidth',5,'TimeWindow',[-1 0.5]);
yline(0,'k--','LineWidth',1.5);xlim([-1 0.5]);xline(0,'k--','LineWidth',1.5);ylim([-0.1 0.1]);
hold off;xlabel('time from choice made');ylabel(['coefficient (ΣV)',newline,'(OT — PL)']);
setFixedPlotArea(ax1, [15 15]);
gcaf1(gca,35);setAxesForPPT(ax1, 36, 2.5);cutBackColor(gca,gcf);
print(gcf,[f_file,'MEG-source-conjunction-interaction-timeCourse.svg'], '-dsvg','-vector','-r600');
%% figure S5D
% interaction
time_idx0=vchan_time>-0.4 & vchan_time<-0.3;
xx2=[mean(smoothdata(squeeze(otsb_chosen1(:,time_idx0,2)),2,'movmean',20),2),...
    mean(smoothdata(squeeze(plsb_chosen1(:,time_idx0,2)),2,'movmean',20),2),...
    mean(smoothdata(squeeze(otnb_chosen1(:,time_idx0,2)),2,'movmean',20),2),...
    mean(smoothdata(squeeze(plnb_chosen1(:,time_idx0,2)),2,'movmean',20),2)];
ax1=figureX1([23 20]);
x_pos=figureTmp1_OT(xx2,1);hold on;plot([xlim],[0 0],'k--','LineWidth',2.5);
setFixedPlotArea(ax1, [15 15]);
box off;ylabel('coefficient (ΣV)');
set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'}); 
gcaf1(gca,35);setAxesForPPT(ax1, 36, 2.5);cutBackColor(gca,gcf);
print(gcf,[f_file,'MEG-source-conjunction-interaction_plot.svg'], '-dsvg','-vector','-r600');



%% MEG source & fMRS overlap
load(fullfile(dataRoot,'MRS_results','MEG_OXTR_fMRS_overlap.mat'),'idx_fmrs','idx_oxtr','idx_meg_tfr','idx_meg_efr');
idx_three_overlap1 = intersect(intersect(idx_meg_tfr,idx_oxtr),idx_fmrs); % meg_gamma-band source ∩ oxtr ∩ fMRS
idx_three_overlap2 = intersect(intersect(idx_meg_efr,idx_oxtr),idx_fmrs); % meg_coefficient(ΣV) source ∩ oxtr ∩ fMRS

%% figure S5E
% S5E left use connectome workbench generate meg_gamma-band source ∩ oxtr ∩ fMRS map
% right: proportion of meg_gamma-band source ∩ oxtr ∩ fMRS
xx1=[length(idx_three_overlap1),length(idx_fmrs)-length(idx_three_overlap1)];
m4_color=[136 73 167; 204 164 225; 221 196 237; 236 223 245]./255;
ax1=figureX1([20 20]);
p_pi=pie(xx1,[1 0]);
for i=1:length(p_pi)./2
    p_pi(2*i-1).FaceColor=m4_color(i,:);
    p_pi(2*i).Position=p_pi(2*i).Position.*0.5;
    p_pi(2*i).Color=[1 1 1];
    p_pi(2*i).FontName='Arial';p_pi(2*i).FontSize=25;
end
setFixedPlotArea(ax1, [15 15]);gcaf1(gca,35);
setAxesForPPT(ax1, 35, 2.5);
print(gcf,[f_file,'conjunction_tfrMEG_fMRS29.svg'], '-dsvg','-vector','-r600');

%% figure S5F
% S5E left use connectome workbench generate meg_coefficient(ΣV) source ∩ oxtr ∩ fMRS map
% right: proportion of meg_coefficient(ΣV) source ∩ oxtr ∩ fMRS
xx2=[length(idx_three_overlap2),length(idx_fmrs)-length(idx_three_overlap2)];
ax1=figureX1([20 20]);
p_pi=pie(xx2,[1 0]);
for i=1:length(p_pi)./2
    p_pi(2*i-1).FaceColor=m4_color(i,:);
    p_pi(2*i).Position=p_pi(2*i).Position.*0.5;
    p_pi(2*i).Color=[1 1 1];
    p_pi(2*i).FontName='Arial';p_pi(2*i).FontSize=25;
end
setFixedPlotArea(ax1, [15 15]);gcaf1(gca,35);
setAxesForPPT(ax1, 35, 2.5);
print(gcf,[f_file,'conjunction_efrMEG_fMRS29.svg'], '-dsvg','-vector','-r600');




%%
%%%%%%%%%%%%%%%%%%%%%%%% function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = bootstrap_trajectory_sem(X,Y,time_idx,varargin)
% X/Y: subject × complete_time
%
% Center trajectory corresponds exactly to: 
% mean across subjects
% → first smooth
% → select time window
% → abs
% → second smooth

p = inputParser;

addParameter(p,'NBoot',5000);
addParameter(p,'FirstSmoothMethod','movmean');
addParameter(p,'FirstSmoothWindow',20);
addParameter(p,'SecondSmoothMethod','gaussian');
addParameter(p,'SecondSmoothWindow',20);
addParameter(p,'Seed',2026);

parse(p,varargin{:});
opt = p.Results;

X = double(X);
Y = double(Y);

assert(isequal(size(X),size(Y)), ...
    'X and Y sizes differ.');

assert(size(X,2)==numel(time_idx), ...
    'time_idx length differs from data time dimension.');

n_sub = size(X,1);
n_selected = nnz(time_idx);

%% Original sample center trajectory
mean_x_signed = mean(X,1,'omitnan');
mean_y_signed = mean(Y,1,'omitnan');

mean_x_signed = smoothdata( ...
    mean_x_signed,2, ...
    opt.FirstSmoothMethod, ...
    opt.FirstSmoothWindow, ...
    'omitnan');

mean_y_signed = smoothdata( ...
    mean_y_signed,2, ...
    opt.FirstSmoothMethod, ...
    opt.FirstSmoothWindow, ...
    'omitnan');

% Identical to the original plotting code: Take absolute values within the selected window, then smooth
mean_x = smoothdata( ...
    abs(mean_x_signed(time_idx)),2, ...
    opt.SecondSmoothMethod, ...
    opt.SecondSmoothWindow, ...
    'omitnan');

mean_y = smoothdata( ...
    abs(mean_y_signed(time_idx)),2, ...
    opt.SecondSmoothMethod, ...
    opt.SecondSmoothWindow, ...
    'omitnan');

%% Paired participants bootstrap
boot_x = nan(opt.NBoot,n_selected,'single');
boot_y = nan(opt.NBoot,n_selected,'single');

rng(opt.Seed,'twister');

for i_boot = 1:opt.NBoot

    % self and other Use exactly the same participant indices
    sub_idx = randi(n_sub,n_sub,1);

    bx = mean(X(sub_idx,:),1,'omitnan');
    by = mean(Y(sub_idx,:),1,'omitnan');

    bx = smoothdata( ...
        bx,2,opt.FirstSmoothMethod, ...
        opt.FirstSmoothWindow,'omitnan');

    by = smoothdata( ...
        by,2,opt.FirstSmoothMethod, ...
        opt.FirstSmoothWindow,'omitnan');

    bx = smoothdata( ...
        abs(bx(time_idx)),2, ...
        opt.SecondSmoothMethod, ...
        opt.SecondSmoothWindow,'omitnan');

    by = smoothdata( ...
        abs(by(time_idx)),2, ...
        opt.SecondSmoothMethod, ...
        opt.SecondSmoothWindow,'omitnan');

    boot_x(i_boot,:) = single(bx);
    boot_y(i_boot,:) = single(by);
end

%% x and y Each one's bootstrap SEM
sem_x = std(boot_x,0,1,'omitnan');
sem_y = std(boot_y,0,1,'omitnan');

%% Local normal to the mean trajectory
% Used only to determine normal directions, Does not change the actual mean trajectory
geometry_n = min(5,n_selected);

geometry_x = smoothdata( ...
    mean_x,2,'movmean',geometry_n,'omitnan');

geometry_y = smoothdata( ...
    mean_y,2,'movmean',geometry_n,'omitnan');

dx = gradient(geometry_x);
dy = gradient(geometry_y);

denom = hypot(dx,dy);

normal_x = nan(1,n_selected);
normal_y = nan(1,n_selected);

valid_normal = isfinite(denom) & denom>eps;

normal_x(valid_normal) = ...
    -dy(valid_normal)./denom(valid_normal);

normal_y(valid_normal) = ...
     dx(valid_normal)./denom(valid_normal);

%% Fill invalid normals
for k = 1:n_selected

    if ~isfinite(normal_x(k)) || ...
            ~isfinite(normal_y(k))

        if k>1
            normal_x(k) = normal_x(k-1);
            normal_y(k) = normal_y(k-1);
        else
            normal_x(k) = 0;
            normal_y(k) = 1;
        end
    end
end

%% Prevent adjacent normals from flipping
for k = 2:n_selected

    if normal_x(k)*normal_x(k-1) + ...
            normal_y(k)*normal_y(k-1)<0

        normal_x(k) = -normal_x(k);
        normal_y(k) = -normal_y(k);
    end
end

%% Bootstrap Project deviations onto the normal
boot_normal = nan(opt.NBoot,n_selected,'single');

for k = 1:n_selected

    boot_normal(:,k) = ...
        (boot_x(:,k)-mean_x(k))*normal_x(k) + ...
        (boot_y(:,k)-mean_y(k))*normal_y(k);
end

% Bootstrap SEM, No further division by sqrt(n_sub)
sem_normal = std(boot_normal,0,1,'omitnan');

%% SEM boundary
upper_x = mean_x + sem_normal.*normal_x;
upper_y = mean_y + sem_normal.*normal_y;

lower_x = mean_x - sem_normal.*normal_x;
lower_y = mean_y - sem_normal.*normal_y;

%% Output
out = struct();

out.mean_x = mean_x;
out.mean_y = mean_y;

out.sem_x = double(sem_x);
out.sem_y = double(sem_y);
out.sem_normal = double(sem_normal);

out.normal_x = normal_x;
out.normal_y = normal_y;

out.upper_x = upper_x;
out.upper_y = upper_y;
out.lower_x = lower_x;
out.lower_y = lower_y;

out.boot_x = boot_x;
out.boot_y = boot_y;

%% SEM Diagnostics relative to the mean trajectory
radius = hypot(mean_x,mean_y);

valid_ratio = ...
    isfinite(radius) & radius>1e-6 & ...
    isfinite(sem_normal);

out.sem_to_radius = nan(size(radius));
out.sem_to_radius(valid_ratio) = ...
    sem_normal(valid_ratio)./radius(valid_ratio);

out.summary = struct();

out.summary.median_sem = ...
    median(sem_normal,'omitnan');

out.summary.max_sem = ...
    max(sem_normal,[],'omitnan');

out.summary.median_radius = ...
    median(radius,'omitnan');

out.summary.median_sem_to_radius = ...
    median(out.sem_to_radius,'omitnan');

out.summary.max_sem_to_radius = ...
    max(out.sem_to_radius,[],'omitnan');

end


function h = draw_trajectory_sem_union(ax,upper_x,upper_y,lower_x,lower_y,sem_color,face_alpha,clip_at_zero)
% DRAW_TRAJECTORY_SEM_UNION Draw a curved 2-D SEM ribbon without alpha buildup.
%
% Each adjacent pair of upper/lower boundary samples defines one local
% quadrilateral. Their geometric union is calculated before plotting, so a
% folded trajectory neither bridges distant parts of the path nor darkens
% where its SEM ribbon overlaps itself.

if nargin<8 || isempty(clip_at_zero)
    clip_at_zero = false;
end

validateattributes(ax,{'matlab.graphics.axis.Axes'},{'scalar'},mfilename,'ax');
upper_x = double(upper_x(:)');
upper_y = double(upper_y(:)');
lower_x = double(lower_x(:)');
lower_y = double(lower_y(:)');
assert(isequal(size(upper_x),size(upper_y),size(lower_x),size(lower_y)), ...
    'All SEM-boundary vectors must have the same size.');
assert(numel(upper_x)>=2,'At least two time samples are required.');
validateattributes(sem_color,{'numeric'}, ...
    {'real','finite','vector','numel',3,'>=',0,'<=',1},mfilename,'sem_color');
validateattributes(face_alpha,{'numeric'}, ...
    {'real','finite','scalar','>=',0,'<=',1},mfilename,'face_alpha');

ribbon_shape = polyshape();
n_time = numel(upper_x);
for k = 1:n_time-1
    vx = [upper_x(k),upper_x(k+1),lower_x(k+1),lower_x(k)];
    vy = [upper_y(k),upper_y(k+1),lower_y(k+1),lower_y(k)];
    if ~all(isfinite([vx,vy]))
        continue;
    end
    segment_shape = polyshape(vx,vy,'Simplify',true);
    if isempty(segment_shape.Vertices) || area(segment_shape)<=eps
        continue;
    end
    if isempty(ribbon_shape.Vertices)
        ribbon_shape = segment_shape;
    else
        ribbon_shape = union(ribbon_shape,segment_shape);
    end
end

if isempty(ribbon_shape.Vertices)
    h = gobjects(0);
    warning('No finite, non-degenerate SEM ribbon could be constructed.');
    return;
end

if clip_at_zero
    vertices = ribbon_shape.Vertices;
    xmax = max(vertices(:,1),[],'omitnan');
    ymax = max(vertices(:,2),[],'omitnan');
    margin = max([xmax,ymax,1])*1e-9;
    clip_shape = polyshape( ...
        [0,xmax+margin,xmax+margin,0], ...
        [0,0,ymax+margin,ymax+margin]);
    ribbon_shape = intersect(ribbon_shape,clip_shape);
end

h = plot(ax,ribbon_shape, ...
    'FaceColor',double(sem_color(:)'), ...
    'FaceAlpha',face_alpha, ...
    'EdgeColor','none', ...
    'HandleVisibility','off');
end

function plot_gradient_path(ax,x,y,timeColors,isDashed,lineWidth)
    n = numel(x);
    if any(~isfinite(x)) || any(~isfinite(y))
        return;
    end
    for k = 1:(n-1)
        if isDashed && mod(floor((k-1)/4),2) == 1
            continue;
        end
        plot(ax,x(k:k+1),y(k:k+1),'-', ...
            'Color',timeColors(k,:),'LineWidth',lineWidth, ...
            'HandleVisibility','off');
    end
    scatter(ax,x(1),y(1),280,'o','MarkerFaceColor',[1 1 1], ...
        'MarkerEdgeColor',timeColors(1,:),'LineWidth',1,'HandleVisibility','off');
    scatter(ax,x(end),y(end),280,'o','MarkerFaceColor',timeColors(end,:), ...
        'MarkerEdgeColor',timeColors(end,:),'LineWidth',1,'HandleVisibility','off');
end


%%



