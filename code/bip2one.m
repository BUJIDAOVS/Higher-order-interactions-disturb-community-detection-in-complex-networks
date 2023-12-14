function [all_in_one] = bip2one(bip)
% BIP2ONE Converts a bipartite network into a unilayer network of one entity type.
%   This function projects a bipartite network (e.g., scientists and papers) onto
%   a unilayer network of one entity type, such as creating a scientists' collaboration
%   network from a scientists-papers bipartite network. The result is a network
%   where nodes represent one type of entity, and edges indicate their collaboration
%   or connection inferred from the bipartite relations.

%% Convert to a Unilayer Network
tic
num_paper = max(bip(:,1)); % Total number of papers
num_author = max(bip(:,2)); % Total number of authors
p_a = sparse(bip(:,1), bip(:,2), ones(1, length(bip)), num_paper, num_author); % Convert list to sparse matrix

auth_coop_nw = p_a' * p_a; % Scientists' collaboration network: Papers co-authored by scientists
auth_coop_nw(1:length(auth_coop_nw) + 1:end) = 0; % Set the main diagonal to zero

paper_cooped_nw = p_a * p_a'; % Papers' collaboration network: Papers co-authored by common scientists
paper_cooped_nw(1:length(paper_cooped_nw) + 1:end) = 0; % Set the main diagonal to zero
disp(['Unilayer network calculation time: ', num2str(toc)]);

%%
all_in_one = auth_coop_nw;

end