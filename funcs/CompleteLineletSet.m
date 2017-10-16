function [ bComplete ] = CompleteLineletSet( linelet_cand, L, conn_map, main_dir )
%COMPLETELINELETSET Summary of this function goes here
%   Detailed explanation goes here

% Linelet connection map
% -- non-zero value of conn_map(i,j) means that linelet i and j are connected such that
%    positive/negative for i is below/above j, and 1/2 for i is left/right to j

bComplete = -1;

% Get the longest linelet
[val_srt, ind_srt] = sort(linelet_cand(:,3), 'descend');

if L.NumObjects == 3  % Compelete
    if sum(conn_map(ind_srt(1),:)~=0) ~= 2, return; end 
    nbs = conn_map(ind_srt(1),:);
    if prod(nbs(nbs~=0)) == -2
        bComplete = 1;        
    end
elseif L.NumObjects > 3
    if abs(val_srt(1) - mean(val_srt(2:end))) <= val_srt(1)/2 
        bComplete = 1;
    else
        bComplete = 0; % Semi-complete
    end
end

% CompleteLineletSet()
% the length of the current linelet set: len
% the number of linelets within the set: nL
% rough slope of current lineletset: mL

% 	the length of the current linelet set: len
% 	the number of linelets within the set: nL
% 	rough slope of current lineletset: mL
% 	if nL >= 3 && diff(mL, mnL) < theta
% 		yes
% 		no
% if nL >= 3 && diff(mL, mnL) < theta
% 	yes
% 	no

end

