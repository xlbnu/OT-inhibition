function setFixedPlotArea(ax, plotSize_cm,rigth_up)
    % Set fixed plot-area dimensions, Use a fixed TightInset value
    % ax: axes handle
    % plotSize_cm: inner plot-area size [width, height](centimeters)
    if nargin<3
    rigth_up=[1.2 0.6];
    end
    
    % Set units to centimeters
    set(ax, 'Units', 'centimeters');
    
    % Use a fixed TightInset value, Rather than obtaining them from the axes
    % Adjust these values as needed
    fixedTightInset = [0.8, 0.6, 0.3, 0.4];  % [Left, Bottom, Right, Top] Units: centimeters
    
    % Compute the total axes size required
    axWidth = plotSize_cm(1) + fixedTightInset(1) + fixedTightInset(3);
    axHeight = plotSize_cm(2) + fixedTightInset(2) + fixedTightInset(4);
    
    % Get figure
    fig = get(ax, 'Parent');
    set(fig, 'Units', 'centimeters');
    
    % Compute axes position
    figPos = get(fig, 'Position');
    figWidth = figPos(3);
    figHeight = figPos(4);
    
    % Align to upper-right corner
    axLeft = figWidth - axWidth - rigth_up(1);  % Align to the right
    axBottom = figHeight - axHeight - rigth_up(2);  % Align to the top
    
    % Set axes position and dimensions
    set(ax, 'Position', [axLeft, axBottom, axWidth, axHeight]);
    
    % Add a window-resize callback
    set(fig, 'SizeChangedFcn', @(~,~) updateFixedPlotArea(ax, plotSize_cm,rigth_up));
end

function updateFixedPlotArea(ax, plotSize_cm,right_up)
    % Update the fixed plot area
    if isvalid(ax)
        % Set units to centimeters
        set(ax, 'Units', 'centimeters');
        
        % Use the same fixed TightInset
        fixedTightInset = [0.8, 0.6, 0.3, 0.4];
        
        % Compute the total axes size required
        axWidth = plotSize_cm(1) + fixedTightInset(1) + fixedTightInset(3);
        axHeight = plotSize_cm(2) + fixedTightInset(2) + fixedTightInset(4);
        
        % Get figure
        fig = get(ax, 'Parent');
        if isvalid(fig)
            set(fig, 'Units', 'centimeters');
            figPos = get(fig, 'Position');
            figWidth = figPos(3);
            figHeight = figPos(4);
            
            % Align to upper-right corner
            axLeft = figWidth - axWidth - rigth_up(1);
            axBottom = figHeight - axHeight - rigth_up(2);
            
            % Set axes position
            set(ax, 'Position', [axLeft, axBottom, axWidth, axHeight]);
        end
    end
end