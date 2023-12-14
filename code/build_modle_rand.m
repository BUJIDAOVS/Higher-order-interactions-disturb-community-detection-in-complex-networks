function [adj_matrix, bip] = build_modle_rand(num_papers, num_scientists, num_community)
% BIP_MODEL_RAND - Builds a randomized bipartite network model
%   This function constructs a randomized adjacency matrix for a bipartite network model representing scientific collaborations.
%   It simulates random collaboration between scientists and papers without specific community structures.
%   Input parameters:
%       num_papers: Number of papers in the network
%       num_scientists: Total number of scientists
%       num_community: Number of scientists in each community (not used in the random model)
%   Output:
%       adj_matrix: Adjacency matrix indicating random collaborations between papers and scientists
%       bip: Sparse matrix format listing non-zero elements of the adjacency matrix

% Define adjacency matrix and scientist IDs for communities A/B
adj_matrix = zeros(num_papers, num_scientists);

% Define probability distribution for the number of collaborators
prob_dist = [0.195836850782712, 0.322045788345430, 0.226145005851651, ...
             0.127770255827746, 0.0637695580140002, 0.0307516835625079, ...
             0.0152161012424174, 0.00928436020895301, 0.00552741617086667, ...
             0.00365297999371530];

% Construct the adjacency matrix
for i = 1:num_papers
    % Randomly select a number of authors for this paper from the probability distribution
    num_authors = randsample(1:10, 1, true, prob_dist);
    authors = randsample(num_scientists, num_authors, false);
    
    % Establish connections between the authors of this paper and other scientists
    adj_matrix(i, authors) = 1;
end

% Find the row and column indices of all non-zero elements
[I, J] = find(adj_matrix);
% Construct a sparse matrix 'bip' listing the links in the network
bip = [I, J];

end
