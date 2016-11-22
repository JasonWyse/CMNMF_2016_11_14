function [ g_p_network ] = Extend_GPNetwork_TruePathRule( g_p_network, M, first_level_id)
%EXTEND_GPNETWORK_TRUEPATHRULE Summary of this function goes here
%   Detailed explanation goes here

    [r,c] = find(M==1);
    for i = 1:length(r)
        row = r(i);
        column = c(i);
        [rows,~] = find(g_p_network(:,length(first_level_id)+column)==1);
        g_p_network(rows,row)=1;
    end

end

