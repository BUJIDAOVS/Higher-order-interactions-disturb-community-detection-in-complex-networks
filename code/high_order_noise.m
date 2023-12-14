%% Model NMI Image
n = 1000;

NMIallorder = zeros(10,n);
NMIallorder2 = zeros(10,n);

for k = 1:n
    % Bipartite network
    % Define the size of the bipartite network
    num_papers = 250;   % Number of papers
    num_scientists = 100;   % Number of scientists
    num_community = 50;  % Number of scientists in each community

    [adj_matrix, bip] = build_modle(num_papers, num_scientists, num_community);

    only_conncomp = 0;
    [A, realS, order] = order2A(bip, num_community, only_conncomp);

    NMI = [];
    for i = 1:length(A)
        sum_mat = zeros(size(A{1}));
        % Loop through the first n cells, summing matrices within the cells
        for j = 1:i
            sum_mat = sum_mat + A{j};
        end
        S = community_louvain(sum_mat); % Cumulative order, weighted
        orderNMI = getNMI(S, realS');
        NMIallorder(order(i+1), k) = orderNMI;

        sum_mat(sum_mat >= 1) = 1; % Unweighted
        S1 = community_louvain(sum_mat); % Cumulative order, unweighted
        orderNMI1 = getNMI(S1, realS');
        NMIallorder1(order(i+1), k) = orderNMI1;

        SS = community_louvain(A{i}); % Non-cumulative order
        orderNMI2 = getNMI(SS, realS');
        NMIallorder2(order(i+1), k) = orderNMI2;
    end
end

% Standard Error Calculation
NMIallorderstd = zeros(10,1);
NMIallorderstd1 = zeros(10,1);
NMIallorderstd2 = zeros(10,1);
NMIallorderse = zeros(10,1);
NMIallorderse1 = zeros(10,1);
NMIallorderse2 = zeros(10,1);

for i = 1:10
    NMIallordermean(i) = sum(NMIallorder(i,:)) / nnz(NMIallorder(i,:));
    NMIallordermean1(i) = sum(NMIallorder1(i,:)) / nnz(NMIallorder1(i,:));
    NMIallordermean2(i) = sum(NMIallorder2(i,:)) / nnz(NMIallorder2(i,:));

    arr = NMIallorder(i,:);
    non_zero_arr = arr(arr ~= 0);  % Extract non-zero elements
    NMIallorderstd(i) = std(non_zero_arr, 0);  % Calculate standard deviation
    NMIallorderse(i) = NMIallorderstd(i) / sqrt(length(non_zero_arr));  % Calculate standard error

    arr1 = NMIallorder1(i,:);
    non_zero_arr1 = arr1(arr1 ~= 0);  % Extract non-zero elements
    NMIallorderstd1(i) = std(non_zero_arr1, 0);  % Calculate standard deviation
    NMIallorderse1(i) = NMIallorderstd1(i) / sqrt(length(non_zero_arr1));  % Calculate standard error

    arr2 = NMIallorder2(i,:);
    non_zero_arr2 = arr2(arr2 ~= 0);  % Extract non-zero elements
    NMIallorderstd2(i) = std(non_zero_arr2, 0);  % Calculate standard deviation
    NMIallorderse2(i) = NMIallorderstd2(i) / sqrt(length(non_zero_arr2));  % Calculate standard error
end

figure
errorbar(NMIallordermean1, NMIallorderse1); % Cumulative order, unweighted
hold on
errorbar(NMIallordermean2, NMIallorderse2); % Non-cumulative order
xlim([1.5,10.5]);
title(sprintf('Model NMI (Number of Scientists: %d, Number of Papers: %d)', num_scientists, num_papers));
xlabel('Order');
ylabel('NMI');
legend('Cumulative Order', 'Non-Cumulative Order');

%%
NMIallordermean = zeros(10,1);
NMIallordermean1 = zeros(10,1);
NMIallordermean2 = zeros(10,1);

NMIallordervar = zeros(10,1);
NMIallordervar1 = zeros(10,1);
NMIallordervar2 = zeros(10,1);

for i = 1:10
    NMIallordermean(i) = sum(NMIallorder(i,:)) / nnz(NMIallorder(i,:));
    NMIallordermean1(i) = sum(NMIallorder1(i,:)) / nnz(NMIallorder1(i,:));
    NMIallordermean2(i) = sum(NMIallorder2(i,:)) / nnz(NMIallorder2(i,:));

    arr = NMIallorder(i,:);
    non_zero_arr = arr(arr ~= 0);  % Extract non-zero elements
    NMIallordervar(i) = var(non_zero_arr);  % Calculate variance
    arr1 = NMIallorder1(i,:);
    non_zero_arr1 = arr1(arr1 ~= 0);  % Extract non-zero elements
    NMIallordervar1(i) = var(non_zero_arr1);  % Calculate variance
    arr2 = NMIallorder2(i,:);
    non_zero_arr2 = arr2(arr2 ~= 0);  % Extract non-zero elements
    NMIallordervar2(i) = var(non_zero_arr2);  % Calculate variance
end

figure
errorbar(NMIallordermean, 3 * NMIallordervar, 'b--'); % Weighted
hold on
errorbar(NMIallordermean1, 3 * NMIallordervar1, 'k--'); % Unweighted
hold on
errorbar(NMIallordermean2, 3 * NMIallordervar2, 'r--');
ylim([0,1]);

%%
%%
NMIallorder = zeros(10, 1000);
for k = 1:1000
    % Bipartite Network
    % Define the size of the bipartite network
    num_papers = 200;   % Number of papers
    num_scientists = 100;   % Number of scientists
    num_community = 50;  % Number of scientists in each community

    [adj_matrix, bip] = build_modle(num_papers, num_scientists, num_community);

    only_conncomp = 0;
    [A, realS, order] = order2A(bip, num_community, only_conncomp);

    NMI = [];
    below_A = cell(1, length(A));
    for i = 1:length(A)
        varname = sprintf('order%d_below', order(i+1));
        sum_mat = zeros(size(A{1}));
        % Loop through the first n cells, adding up the matrices in the cells
        for j = 1:i
            sum_mat = sum_mat + A{i};
        end
        sum_mat(sum_mat >= 1) = 1; % Unweighted
        eval([varname, ' = sum_mat;']);
        below_A{i} = sum_mat;
        S = community_louvain(sum_mat); % Community detection
        orderNMI = getNMI(S, realS'); % Calculate NMI
        NMIallorder(order(i+1), k) = orderNMI;
    end
end

% Calculate the mean NMI for each order
NMIallordermean = zeros(10,1);
for i = 1:10
    NMIallordermean(i) = sum(NMIallorder(i,:)) / nnz(NMIallorder(i,:));
end

% Plot the mean NMI
plot(NMIallordermean);
%%






