function [ID_out] = bot_para_o_da(ID0,ID1)
% Generate a categorical array
% id É O VECTOR COM AS 3 EWPLICAS
% Generate a categorical array


% Set the number of bootstrap samples and sample size
nboot = 200;
n = min(numel(ID0),numel(ID1));

% % Initialize matrix for bootstrap samples
% ID_0 = categorical(zeros(n, nboot));

% Generate bootstrap samples
for i = 1:nboot
    ID_0(:,1) = datasample(ID0, n, 'Replace', true);
    ID_1(:,1)=  datasample(ID1, n, 'Replace', true);
    ID_out(:,i) =[ID_0;ID_1];
end

end