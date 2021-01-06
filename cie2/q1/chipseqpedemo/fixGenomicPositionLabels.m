function fixGenomicPositionLabels(ha)
% FIXGENOMICPOSITIONLABELS Helper function that adds callback functions for
% updating XTickLabel and datacursor, so the X axis labels are formatted
% properly as genomic positions. FIXGENOMICPOSITIONLABELS also turns the
% box on and resizes the figure, so the rendering is more appealing to
% genomic data.
%
% FIXGENOMICPOSITIONLABELS(HA) HA is the axes containing the genomic data.

%   Copyright 2012-2015 The MathWorks, Inc.

if nargin == 0
    ha = gca;
end

hf = get(ha,'Parent');   % handle to figure
dc = datacursormode(hf); % handle to datacursor manager

ad.UpdateFunctionToHonor = get(dc,'UpdateFcn'); 
setappdata(hf,'fixGenomicPositionLabels',ad)  % store old datacursor function handle
set(dc,'UpdateFcn',@updateDataCursor) 

datacursormode(hf,'on')
setappdata(ha,'datacursorHandledByFixGenomicPositionLabels',true);

lh = addlistener(ha, 'SizeChanged',@updateXTickLabels);
setappdata(ha,'XtickListener',lh);

% turn on box and resize the window
box(ha,'on')             
set(hf,'Position',max(get(hf,'Position'),[0 0 900 0])) 

end

function updateXTickLabels(~,e)
% Callback function to reformat the XTick labels

ha = e.Source;
set(ha,'XTickLabelMode','auto');

% format the labels 
v = get(ha,'XTick');
labels = arrayfun(@(x)(num2str(x)),v,'UniformOutput',false);
labels = addCommasToLabels(labels);

% find out which labels are not displayed in auto mode and set them to ''
oldlabels = get(ha,'xticklabel');
idx = cellfun(@(x)(isempty(x)),oldlabels); 
labels(idx) = {''}; 

% set the newly formatted labels
set(ha,'XTickLabel',labels);

end

function txt = updateDataCursor(h,e)
% Callback function to prepare the text that goes in a datacursor

% find the first parent axes to the target (calling target may be inside hggropus)
ha = get(e,'Target');
while ~strcmp(get(ha,'Type'),'axes') && ha~=0
    ha = get(ha,'Parent');
end

% Check if fixGenomicPositionLabels owns the data cursor for this axes
if ~isempty(getappdata(ha,'datacursorHandledByFixGenomicPositionLabels'))

    % prepare the datacursor for this axes
      posLabelsX = addCommasToLabels({num2str(get(e,'Position')*[1;0])});
      posLabelsY = num2str(get(e,'Position')*[0;1]);
      txt = {['Position: ',posLabelsX{:}],['Y: ', posLabelsY]};
else
    hf = get(ha,'Parent'); %handle to figure
    ad = getappdata(hf,'fixGenomicPositionLabels');
    
    % check if there is another updateDataCursor that neeeds to be called
    if ~isempty(ad.UpdateFunctionToHonor)
        % Honor the previously saved UpdateFcn in axes that I do not own
        txt = feval(ad.UpdateFunctionToHonor,h,e);
    else
        % fixGenomicPositionLabels does not own the axes and the
        % previously saved UpdateFcn is empty, then just replicate MATLAB's
        % default datacursor
        pos = get(e,'position');
        txt = {['X=',num2str(pos(1))],['Y=',num2str(pos(2))]};
    end
end
end

function fLabels = addCommasToLabels(labels)
% Format strings so that commas are inserted every 3 digits from left-most
% positions (e.g. 2000000 -> 2,000,000)

fLabels = cellfun(@(x)(fliplr(x)),labels,'UniformOutput',false);
fLabels = cellfun(@(x)(sprintf(',%c%c%c',x)),fLabels,'UniformOutput',false);
fLabels = cellfun(@(x)(fliplr(x(2:end))),fLabels,'UniformOutput',false);

end