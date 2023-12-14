% Bipartite network to normal network
% Convert paper-author bipartite network to a unilayer scientist collaboration network

%% Clear workspace
clc; clear;

%% Set working directory
cd("/Users/liuyuyan/Documents/MATLAB/projects/sci_coop_nw/");

% Uncomment to read data from nw_a (Large network APS American Physical Science magazine)
% ...

%% Reading data from nw_b (Small network Econophysics)
bip = readmatrix("./data/nw_b/Scientists_collaboration_networks/paperauthorId");
bip = bip + 1; % Adjusting IDs to start from 1 for MATLAB indexing

%% Convert to a Unilayer Network
tic
num_paper = max(bip(:,1)); % Total number of papers
num_author = max(bip(:,2)); % Total number of authors
p_a = sparse(bip(:,1), bip(:,2), ones(1, length(bip)), num_paper, num_author); % Convert list to sparse matrix

auth_coop_nw = p_a' * p_a; % Scientist collaboration network: Scientists co-authoring papers
auth_coop_nw(1:length(auth_coop_nw)+1:end) = 0; % Zero out the main diagonal

paper_cooped_nw = p_a * p_a'; % Paper collaboration network: Papers co-authored by common scientists
paper_cooped_nw(1:length(paper_cooped_nw)+1:end) = 0; % Zero out the main diagonal
disp(['Unilayer network calculation time: ', num2str(toc)]);

%% High-Order Network
% Sorting papers by the number of collaborating authors
frequence = tabulate(bip(:,1));
frequence = sortrows(frequence, 2, 'descend');
order = unique(frequence(:,2));

% Looping through each order
for i = order'
    tic
    % Extracting bipartite network for papers co-authored by 'i' authors
    clear i_order_paper;
    i_order_paper = frequence(frequence(:,2) == i, 1);
    i_order_bip = [];
    for j = i_order_paper'
        i_order_bip = [i_order_bip; bip(bip(:,1) == j, :)];
    end

    % Converting extracted network to i-th order network
    i_order_p_a = sparse(i_order_bip(:,1), i_order_bip(:,2), ones(1, length(i_order_bip)), num_paper, num_author);
    i_order_auth_coop_nw = i_order_p_a' * i_order_p_a;
    i_order_auth_coop_nw(1:length(i_order_auth_coop_nw)+1:end) = 0;

    % Store the i-th order network in a variable named 'orderi_auth_coop_nw'
    i_order_auth_coop_nw_name = ['order', mat2str(i), '_auth_coop_nw'];
    eval([i_order_auth_coop_nw_name, '=i_order_auth_coop_nw;']);

    % Optionally, save the i-th order network
    % eval(['save(','''./data/nw_b/save_eye/', i_order_auth_coop_nw_name, '.mat'',''', i_order_auth_coop_nw_name, ''');']);
    disp([[mat2str(i), ' order network calculation time: '], num2str(toc)]);
end
