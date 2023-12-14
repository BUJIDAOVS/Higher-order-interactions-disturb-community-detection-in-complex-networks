%%

% Reading data from a CSV file
paper_year_pacs = csvread('/Users/liuyuyan/Documents/MATLAB/projects/sci_coop_nw/data/nw_a/paperid_year_pacs_xx.csv');

% Adjusting paper IDs to start from 1 (as MATLAB arrays start from 1)
paper_year_pacs(:,1) = paper_year_pacs(:,1) + 1;

% Extracting paper IDs and years
paper = paper_year_pacs(:,1);
year = paper_year_pacs(:,2);

% Extracting PACS codes
pacs = paper_year_pacs(:,3:end);

% Reading author-paper relationship data
authorid_paperid = csvread('/Users/liuyuyan/Documents/MATLAB/projects/sci_coop_nw/data/nw_a/authorid_paperid.csv');
bip = readmatrix("/Users/liuyuyan/Documents/MATLAB/projects/sci_coop_nw/data/nw_a/authorid_paperid.csv");

% Adjusting IDs to start from 1 for MATLAB compatibility
bip = bip + 1; 

% Analyzing paper count by year
paper_by_year = [];
for i = 1893:2010
    paper_by_year(i - 1892) = length(intersect(paper(find(year == i)), bip(:,2)));
end

% Plotting the number of papers by year
figure
plot([1893:2010], paper_by_year);

%%
year2 = [];
for i = 1:length(bip(:,2))
    if find(paper == bip(i,2)) ~= 0
        year2(i) = year(find(paper == bip(i,2)));
    else 
        year2(i) = nan;
    end
end

save('C:\Users\bujidao\Documents\MATLAB\paper_year\year2.mat',"year2")

%% Preprocessing - Formatting data as authorid_paperid_year_pacs

year_pacs= [];
for i = 1:length(bip(:,2))
    if find(paper == bip(i,2)) ~= 0
        year_pacs(i,:) = paper_year_pacs(find(paper == bip(i,2)),:);
    else
        year_pacs(i,:) = nan;
    end
end

authorid_paperid_year_pacs = [bip,year_pacs(:,2:end)];

%%
save('/Users/liuyuyan/Documents/MATLAB/projects/sci_coop_nw/data/nw_a/authorid_paperid_year_pacs.mat',"authorid_paperid_year_pacs")

%%
csvwrite('/Users/liuyuyan/Documents/MATLAB/projects/sci_coop_nw/data/nw_a/authorid_paperid_year_pacs.csv',authorid_paperid_year_pacs)

%%
% Loading data from the specified file
load('/Users/liuyuyan/Documents/MATLAB/projects/sci_coop_nw/data/nw_a/authorid_paperid_year_pacs.mat')

% Assigning loaded data to a variable
apyp = authorid_paperid_year_pacs;
% Converting NaN values to 0
apyp(find(isnan(apyp) == 1)) = 0;

% Retaining valid data
% Extracting rows where the third and fourth columns are non-zero
apyp = apyp(apyp(:,3) ~= 0 & apyp(:,4) ~= 0, :);

% Removing papers published before 1999
apyp = apyp((apyp(:,3) >= 1999), :);

% Creating a bipartite network
bip = apyp(:,1:2);
tem = bip(:,2);
bip(:,2) = bip(:,1);
bip(:,1) = tem; % Swapping columns to format from author-paper to paper-author
clear tem;

% Condensing into a single-layer network
all_in_one = bip2one(bip);

% Finding the largest connected component and saving its index
[com_pacs, sizes] = conncomp(graph(all_in_one));
% Identifying the largest connected component
[~, idx] = max(sizes);

% Initializing a sparse matrix L
L = sparse(size(all_in_one));

% Assigning the largest connected component to L
L(com_pacs == idx, com_pacs == idx) = all_in_one(com_pacs == idx, com_pacs == idx);

%%
S = community_louvain(L);


%% Community detection in the weighted network L
gamma = 1;
k = full(sum(L));
twom = sum(k);
B = @(i) L(:,i) - gamma*k'*k(i)/twom;
[S,Q] = genlouvain(B);
Q = Q/twom;

%% Community detection in the unweighted network L
L_to_one = L;
L_to_one(L_to_one ~= 0) = 1;
gamma = 1;
k = full(sum(L_to_one));
twom = sum(k);
B = @(i) L_to_one(:,i) - gamma*k'*k(i)/twom;
[S_to_one,Q] = genlouvain(B);
Q = Q/twom;
%% 
%% Conversion of bipartite network into hierarchical networks
% Sorting papers based on the number of co-authors
num_paper = max(bip(:,1)); % Total number of papers
num_author = max(bip(:,2)); % Total number of authors
frequence = tabulate(bip(:,1));
frequence = sortrows(frequence,2,'descend');
order = unique(frequence(:,2));
for i = order'
    tic
    % Extracting the bipartite network of papers co-authored by 'i' authors (the network of author collaborations on 'i-author papers' is referred to as the 'i-th order network')
    clear i_order_paper;
    i_order_paper = frequence(frequence(:,2) == i,1);
    i_order_bip = [];
    for j = i_order_paper'
        i_order_bip = [[i_order_bip];[bip(bip(:,1)==j,:)]];
    end

    % Converting the extracted bipartite network of papers co-authored by 'i' authors into the i-th order network
    i_order_p_a = sparse(i_order_bip(:,1),i_order_bip(:,2),ones(1,length(i_order_bip)),num_paper,num_author); % Convert list to sparse matrix
    i_order_auth_coop_nw = i_order_p_a' * i_order_p_a; % Scientist collaboration network where scientists co-author papers
    i_order_auth_coop_nw(1:length(i_order_auth_coop_nw)+1:end) = 0; % Zeroing the main diagonal
    
    % Storing the obtained i-th order network in a variable named 'orderi_auth_coop_nw'
    i_order_auth_coop_nw_name = ['order',mat2str(i),'_auth_coop_nw'];
    eval([i_order_auth_coop_nw_name,'=i_order_auth_coop_nw;']);
    
    % Uncomment to save the network: 
    % eval(['save(','''./data/nw_b/save_eye/',i_order_auth_coop_nw_name,'.mat'',''',i_order_auth_coop_nw_name,''');']); % Save
    disp([[mat2str(i),'th order network computation time: '],num2str(toc)]);
end


%% Association between community IDs and PACS codes
com_num = size(unique(S),1);
com_pacs = cell(com_num, 1);
for i = 1:length(com_pacs)
    com_pacs{i} = [];
end

% authorid_community = [[1:size(S,1)]',S];

for i = [1:size(S,1)]
    tem = apyp(find(apyp(:,1)==i),4:end);
    com_pacs{S(i)} = [com_pacs{S(i)},tem(:)'];
    com_pacs{S(i)}(com_pacs{S(i)}==0)=[];
    clear tem;
end


%% Association between community IDs and PACS codes (all orders, excluding minor authors)
com_num = size(unique(S),1);
com_pacs = cell(com_num, 1);
for i = 1:length(com_pacs)
    com_pacs{i} = [];
end

% authorid_community = [[1:size(S,1)]',S];

for i = [1:size(S,1)]
    if length(find(apyp(:,1)==i)) >= 30
        tem = apyp(find(apyp(:,1)==i),4:end);
        com_pacs{S(i)} = [com_pacs{S(i)},tem(:)'];
        com_pacs{S(i)}(com_pacs{S(i)}==0)=[];
        clear tem;
    end
end

%%


order_more_then_15 = order0_auth_coop_nw;
for i = order'
    if i > 15
        i_order_auth_coop_nw_name = ['order',mat2str(i),'_auth_coop_nw'];
        eval(['order_more_then_15 = order_more_then_15+',i_order_auth_coop_nw_name,';']);
    end
end

all_in_one_less_then_15 = all_in_one - order_more_then_15;

%%
L_less_then_15 = sparse(size(all_in_one_less_then_15));
L_less_then_15(com_pacs == idx, com_pacs == idx) = all_in_one_less_then_15(com_pacs == idx, com_pacs == idx);


%% Network of orders less than or equal to 15 and its largest connected component
order_less_then_15 = order0_auth_coop_nw;
for i = order'
    if i <= 15
        i_order_auth_coop_nw_name = ['order',mat2str(i),'_auth_coop_nw'];
        eval(['order_less_then_15 = order_less_then_15+',i_order_auth_coop_nw_name,';']);
    end
end

% Consolidating the network for orders less than or equal to 15
all_in_one_less_then_15 = order_less_then_15;

% Finding the largest connected component and saving its index
[com_pacs, sizes] = conncomp(graph(all_in_one_less_then_15));
% Identifying the largest connected component
[~, idx] = max(sizes);

% Initializing a sparse matrix for the largest connected component
L_less_then_15 = sparse(size(all_in_one_less_then_15));
% Assigning the largest connected component to L_less_then_15
L_less_then_15(com_pacs == idx, com_pacs == idx) = all_in_one_less_then_15(com_pacs == idx, com_pacs == idx);


%% Community detection in networks of orders less than or equal to 15

gamma = 1;
k = full(sum(L_less_then_15));
twom = sum(k);
B_less_then_15 = @(i) L_less_then_15(:,i) - gamma*k'*k(i)/twom;
[S_less_then_15,Q] = genlouvain(B_less_then_15);
Q = Q/twom;

%% Association between community IDs and PACS codes for networks less than 15
com_num_less_then_15 = size(unique(S_less_then_15),1);
com_pacs_less_then_15 = cell(com_num_less_then_15, 1);
for i = 1:length(com_pacs_less_then_15)
    com_pacs_less_then_15{i} = [];
end

for i = [1:size(S_less_then_15,1)]
    tem = apyp(find(apyp(:,1)==i),4:end);
    com_pacs_less_then_15{S_less_then_15(i)} = [com_pacs_less_then_15{S_less_then_15(i)},tem(:)'];
    com_pacs_less_then_15{S_less_then_15(i)}(com_pacs_less_then_15{S_less_then_15(i)}==0)=[];
    clear tem;
end

%% Association between community IDs and PACS codes for networks less than 15 (excluding minor authors)
com_num_less_then_15 = size(unique(S_less_then_15),1);
com_pacs_less_then_15 = cell(com_num_less_then_15, 1);
for i = 1:length(com_pacs_less_then_15)
    com_pacs_less_then_15{i} = [];
end

for i = [1:size(S_less_then_15,1)]
    if length(find(apyp(:,1)==i)) >= 30
        tem = apyp(find(apyp(:,1)==i),4:end);
        com_pacs_less_then_15{S_less_then_15(i)} = [com_pacs_less_then_15{S_less_then_15(i)},tem(:)'];
        com_pacs_less_then_15{S_less_then_15(i)}(com_pacs_less_then_15{S_less_then_15(i)}==0)=[];
        clear tem;
    end
end

%%

% Creating a cell array
% C = com_pacs_less_then_15;
C = com_pacs;

% Using the cellfun function with the isempty function to get the logical values of non-empty array elements
notEmpty = cellfun(@(x) ~isempty(x), C);

% Calculating the number of non-empty array elements
count = sum(notEmpty)

%%
% Creating a cell array
% C = com_pacs_less_then_15;
C = com_pacs;

% Using the cellfun function to get logical indices of non-empty elements
idx = cellfun(@isempty, C);
idx = ~idx;

% Displaying the indices of non-empty elements
disp(find(idx));


%% Community authors' field similarity (calculated using networks of orders less than or equal to 15)
clear sim;

% Creating a cell array
C = com_pacs_less_then_15;
% C = com_pacs;

% Using cellfun function to get logical indices of non-empty elements
idx = cellfun(@isempty, C);
idx = ~idx;

% Calculating frequency distribution of community IDs in the network
fre = tabulate(S_less_then_15);
% Selecting community IDs with more than one member
com_id_useful = fre((find(fre(:,2) > 1)),1);
% com_id_useful = find(idx);

% Calculating similarity for each selected community
for j = 1:length(com_id_useful)
    clear tem_fre i;
    i = com_id_useful(j);
    if isempty(com_pacs_less_then_15{i})
        continue
    end
    % Tabulating PACS codes in the community
    tem_fre = tabulate(fix(com_pacs_less_then_15{i}));
    % Sorting frequencies in descending order
    sorted = sort(tem_fre(:,3), 'descend');
    % Calculating similarity as the proportion of the most frequent PACS code
    sim(i) = sum(sorted(1:1))/100;
end
% Removing zeros from the similarity array
sim(sim == 0) = [];
% Calculating the mean similarity across all communities
mean(sim)

%% Community authors' field similarity (calculated using networks of all orders)
clear sim;

% Creating a cell array
% C = com_pacs_less_then_15;
C = com_pacs;

% Using cellfun function to get logical indices of non-empty elements
idx = cellfun(@isempty, C);
idx = ~idx;

% Calculating frequency distribution of community IDs in the network
fre = tabulate(S);
% Selecting community IDs with more than one member and excluding all single-member communities
com_id_useful = fre((find(fre(:,2) > 1)),1);
% com_id_useful = find(idx);   % Retaining single-member communities with PACS codes

% Calculating similarity for each selected community
for j = 1:length(com_id_useful)
    clear tem_fre i;
    i = com_id_useful(j);
    if isempty(com_pacs{i})
        continue
    end
    % Tabulating PACS codes in the community
    tem_fre = tabulate(fix(com_pacs{i}));
    % Sorting frequencies in descending order
    sorted = sort(tem_fre(:,3), 'descend');
    % Calculating similarity as the proportion of the most frequent PACS code
    sim(i) = sum(sorted(1:1))/100;
end
% Removing zeros from the similarity array
sim(sim == 0) = [];
% Calculating the mean similarity across all communities
mean(sim)



%% Number of members within each community
for j = 1:length(com_id_useful)
    clear tem_fre i;
    i = com_id_useful(j);
    comid_authnum(j,:) = [i,length(find(S == i))];
    
end
