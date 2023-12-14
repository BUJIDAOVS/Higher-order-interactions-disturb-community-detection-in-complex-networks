function [adj_matrix, bip] = build_modle(num_papers, num_scientists, num_community)
% BIP_MODEL - Builds a bipartite network model
%   This function constructs an adjacency matrix for a bipartite network model representing scientific collaborations.
%   It simulates the collaboration between scientists and papers within or across different communities.
%   Input parameters:
%       num_papers: Number of papers in the network
%       num_scientists: Total number of scientists
%       num_community: Number of scientists in each community
%   Output:
%       adj_matrix: Adjacency matrix indicating collaborations between papers and scientists
%       bip: A sparse matrix format listing non-zero elements of the adjacency matrix, representing the links in the network.
%
%   The function uses a probability distribution to randomly assign authors to papers, considering both within-community and across-community collaborations.


% Define adjacency matrix and scientist IDs for communities A/B
adj_matrix = zeros(num_papers, num_scientists);
community_a = 1:num_community;
community_b = num_community+1:num_scientists;

% Define probability distribution for the number of collaborators
prob_dist = [0.195836850782712, 0.322045788345430, 0.226145005851651, ...
             0.127770255827746, 0.0637695580140002, 0.0307516835625079, ...
             0.0152161012424174, 0.00928436020895301, 0.00552741617086667, ...
             0.00365297999371530];

team_inside_paperid = [];
team_across_paperid = [];

% Construct the adjacency matrix
for i = 1:num_papers
    % Randomly choose a number from the probability distribution to determine the number of authors for this paper
    num_authors = randsample(1:10, 1, true, prob_dist);

    % Probability criterion for determining if the paper is a within-community or across-community collaboration
    prob_line = num_authors/15 - 1/15;
    if rand() > prob_line
        % Collaboration within one community
        team_inside_paperid = [team_inside_paperid, i];
        if rand() < 0.5
            % Select num_authors scientists from community A as authors of the paper
            authors = randsample(community_a, num_authors, false);
        else
            % Select num_authors scientists from community B as authors of the paper
            authors = randsample(community_b, num_authors, false);
        end
    else
        % Collaboration across communities
        team_across_paperid = [team_across_paperid, i];
        if mod(num_authors, 2) == 0
            a_authors = randsample(community_a, floor(num_authors/2), false);
        else
            if rand() < 0.5
                a_authors = randsample(community_a, floor(num_authors/2), false);
            else
                a_authors = randsample(community_a, floor(num_authors/2) + 1, false);
            end
        end
        b_authors = randsample(community_b, num_authors - length(a_authors), false);
        authors = [a_authors, b_authors];
    end
    
    % Establish connections between the authors of this paper and other scientists
    adj_matrix(i, authors) = 1;
end

[I, J] = find(adj_matrix); % Find the row and column indices of all non-zero elements
bip = [I, J]; % Construct a sparse matrix 'bip'

end
