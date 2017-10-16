function [ handles ] = bring_figure_front( hFig  )
%BRING_ALL_FIGURES_FRONT Summary of this function goes here
%   Detailed explanation goes here

if nargin == 0
    handles=findall(0,'type','figure');
    for i = 1:numel(handles)
        figure(handles(i));
    end
else
    figure(hFig);
end

end

