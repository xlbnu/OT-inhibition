function hh = errorbound(varargin)
% ERRORBOUND errorbar like plot, using bound area instead of errorbars.
% Usage:
%       ERRORBOUND(y, e)
%       ERRORBOUND(x, y, e)
%       ERRORBOUND(x, y, l, u)
%       ERRORBOUND(..., 'LineSpec')
%       ERRORBOUND(..., Name, Value)
%
% See also: errorbar
%
% This function is released as part of the `FSLboost` package.
%                                   Powered by Matlab && FSL.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2014-2016 LiTuX @BNU, all rights reserved.
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
[args, pvpairs] = parseparams(varargin);
nargs = length(args);
switch nargs
    case 2
        [y, e] = deal(args{:});
        x = (1: length(y))';
        u = abs(e);
        l = u;
    case 3
        [x, y, e] = deal(args{:});
        u = abs(e);
        l = u;
    case 4
        [x, y, l, u] = deal(args{:});
    otherwise
        error('Error calling %s.', mfilename);
end

colors = 'brgcmyk';
% persistent coloridx;
% if isempty(coloridx)
%     coloridx = 1;
% else
%     coloridx = mod(coloridx, length(colors)) + 1;
% end
option.color = [];      % auto-color
option.linestyle = '-';
option.linewidth = 2;
option.marker = '';
option.facealpha = 0.15;
option.edgestyle = '--';
option.edgewidth = 0.5;

lpropexclude = {'facecolor'; 'facealpha'; 'edgecolor'};
lprops = {};
if ~isempty(pvpairs)
    [linestyle,color,marker,tmsg] = colstyle(pvpairs{1}, 'plot');
    if isempty(tmsg)
        pvpairs = pvpairs(2:end);
        if ~isempty(linestyle)
            pvpairs = [{'LineStyle', linestyle}, pvpairs];
        end
        if ~isempty(color)
            pvpairs = [{'Color',color},pvpairs];
        end
        if ~isempty(marker)
            pvpairs = [{'Marker',marker},pvpairs];
        end
    end
    
    for ii = 1: 2: length(pvpairs)
        prop = lower(pvpairs{ii});
        value = pvpairs{ii+1};
        option.(prop) = value;
        if ~ismember(prop, lpropexclude)
            lprops = cat(2, lprops, {prop, value});
        end
    end
    
end

if min(size(x)) == 1
    x = x(:);
end

if min(size(y)) == 1
    y = y(:);
end
if min(size(u)) == 1
    u = u(:);
end
if min(size(l)) == 1
    l = l(:);
end

if size(x, 1) == 1
    x = repmat(x, 1, size(y,2));
end
if ~isequal(size(x), size(y), size(u), size(l))
    error('Input size miss match.');
end

h = gcf;
holdstatus = ishold;
hold on;
%%
for ii = 1: size(y, 2)
    nanidx = isnan(x(:, ii)) | isnan(y(:, ii));
    segidx = [0; find(nanidx); size(y,1)+1];
    
    if ~isempty(option.color)
        thiscolor = option.color;
    else
        thiscolor = colors( mod(ii-1, length(colors))+1 );
    end
    
    if isfield(option, 'facecolor')
        facecolor = option.facecolor;
    else % auto
        facecolor = thiscolor;
    end
    
    autoec = false;
    if isfield(option, 'edgecolor')
        edgecolor = option.edgecolor;
    else
        autoec = true;
        edgecolor = thiscolor;
    end
    
    % draw the patches and put them in the background
    for jj = 1: length(segidx)-1
        xx = x(segidx(jj)+1: segidx(jj+1)-1, ii);
        yy = y(segidx(jj)+1: segidx(jj+1)-1, ii);
        uu = u(segidx(jj)+1: segidx(jj+1)-1, ii);
        ll = l(segidx(jj)+1: segidx(jj+1)-1, ii);
        
        hx = patch('XData', [xx; flipud(xx)], 'YData', [yy-ll; flipud(yy+uu)], ...
            'FaceColor', facecolor, 'EdgeColor', 'None', 'FaceAlpha', option.facealpha);
        set(get(get(hx, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
        % TODO, inversed order:
        uistack(hx, 'bottom');
    end
    
    % then draw boundaries
    xx = x(:, ii);
    yy = y(:, ii);
    uu = u(:, ii);
    ll = l(:, ii);
    
    if ~strcmpi(edgecolor, 'none')
        hx = plot(xx, [yy-ll, yy+uu], option.edgestyle, 'Color', edgecolor, ...
            'LineWidth', option.edgewidth);
        set(get(get(hx(1), 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
        set(get(get(hx(2), 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
        if autoec
            cc = get(hx(1), 'Color');
            set(hx(1), 'Color', cc+0.55*(1-cc));
            cc = get(hx(2), 'Color');
            set(hx(2), 'Color', cc+0.55*(1-cc));
        end
    end
    
    % and the main line
    plot(xx, yy, option.linestyle, 'Color', thiscolor, ...
        'LineWidth', option.linewidth, lprops{:});
end

%%
if ~holdstatus
    hold off;
end

if nargout == 1
    hh = h;
end
end % function

% End Of File
