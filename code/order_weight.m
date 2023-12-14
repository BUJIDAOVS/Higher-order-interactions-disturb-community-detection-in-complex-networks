%% APS Data Processing
cd("/Users/liuyuyan/Documents/MATLAB/projects/sci_coop_nw/");
bip = readmatrix("./data/nw_a/authorid_paperid.csv");
bip = bip + 1; % Adjusting IDs to start from 1 to match MATLAB's indexing
tem = bip(:,2);
bip(:,2) = bip(:,1);
bip(:,1) = tem; % Swapping columns: paper-author -> author-paper
clear tem;

load('/Users/liuyuyan/Documents/MATLAB/projects/sci_coop_nw/data/nw_a/clear_eye/all_order.mat');

% Frequency calculation and sorting
frequence = tabulate(bip(:,1));
frequence = sortrows(frequence, 2, 'descend');
order = unique(frequence(:,2));

[all_in_one] = bip2one(bip); % Weighted all_in_one

% Connected components calculation
[C, sizes] = conncomp(graph(all_in_one));
[~, idx] = max(sizes);
L = all_in_one(C == idx, C == idx);

% Loop through orders
for i = 2:length(order)
    % Generate variable name based on i
    varname = sprintf('order%d_auth_coop_nw', order(i));
    % Convert variable name to a matrix and store in A{i-1}
    temp = eval(varname);
    A{i-1} = temp; % Isolating nodes
    clear temp;
end

% Plotting
figure
for j = 1:2:5 % Loop through specific orders
    weight = [];
    for i = 1:length(A)
        weight(i) = length(find(all_in_one(find(A{i}>0))>j)) / length(find(A{i}>0));
    end
    plot(order(2:end), weight);
    hold on;
end
xlim([2,60])
xlabel('Order') 
ylabel('Proportion of Edge Weights Greater Than 1 per Order') 
title('APS Top 10')

%% Plotting on a single graph with weights greater than 1

orderx = [1:10];

% Original model to be modified
clear all;
clc;
n = 1000;
weight = zeros(10,n);

for k = 1:n

    num_papers = 250;   % Number of papers
    num_scientists = 100;   % Number of scientists
    num_community = 50;  % Number of scientists in each community

    [adj_matrix,bip] = build_modle(num_papers,num_scientists,num_community);

    only_conncomp = 1;
    [A,realS,order] = order2A_w(bip,num_community,only_conncomp);
    [all_in_one] = bip2one(bip);

    [C, sizes] = conncomp(graph(all_in_one));
    [~, idx] = max(sizes);
    L = all_in_one(C == idx, C == idx);

    for j = 1:1
        for i = 1:length(A)
            weight(order(i+1),k) = length(find(L(find(A{i}>0))>j))/length(find(A{i}>0));
        end
    end
end

% Calculating standard deviation and standard error
weightstd = zeros(10,1);
weightsem = zeros(10,1);

for i = 1:10
    weightmean(i) = sum(weight(i,:))/nnz(weight(i,:));

    arr = weight(i,:);
    non_zero_arr = arr(arr ~= 0);
    weightstd(i) = std(non_zero_arr,0);

    weightsem(i) = weightstd(i) / sqrt(nnz(weight(i,:)));
end

figure
yyaxis left;
orderx = [2:10];
weightmean(1) = [];
weightmean = weightmean*100;
scatter(orderx, weightmean);
hold on;

% Linear fitting
coefficients = polyfit(orderx, weightmean, 1);
fitted_y = polyval(coefficients, orderx);
plot(orderx, fitted_y, 'r');

% Calculating and displaying the slope of the line
slope = coefficients(1);
text(orderx(end), fitted_y(end), sprintf('Slope: %.2f', slope), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

xlim([1.5,10.5]);
ytickformat('percentage')

% New model to be modified
keepvar = {'n'};
whoslist = evalin('base','who');
remvars = setdiff(whoslist, keepvar);
evalin('base', ['clear ' sprintf('%s ', remvars{:})]);

clc;
weight = zeros(10,n);

for k = 1:n

    num_papers = 250;
    num_scientists = 100;
    num_community = 50;

    [adj_matrix,bip] = build_modle_rand(num_papers,num_scientists,num_community);

    only_conncomp = 1;
    [A,realS,order] = order2A_w(bip,num_community,only_conncomp);
    [all_in_one] = bip2one(bip);

    [C, sizes] = conncomp(graph(all_in_one));
    [~, idx] = max(sizes);
    L = all_in_one(C == idx, C == idx);

    for j = 1:1
        for i = 1:length(A)
            weight(order(i+1),k) = length(find(L(find(A{i}>0))>j))/length(find(A{i}>0));
        end
    end
end

% Calculating standard deviation and standard error
weightstd = zeros(10,1);
weightsem = zeros(10,1);

for i = 1:10
    weightmean(i) = sum(weight(i,:))/nnz(weight(i,:));

    arr = weight(i,:);
    non_zero_arr = arr(arr ~= 0);
    weightstd(i) = std(non_zero_arr,0);

    weightsem(i) = weightstd(i) / sqrt(nnz(weight(i,:)));
end

yyaxis left;
hold on
orderx = [2:10];
weightmean(1) = [];
weightmean = weightmean*100;
scatter(orderx, weightmean);
hold on;

% Linear fitting
coefficients = polyfit(orderx, weightmean, 1);
fitted_y = polyval(coefficients, orderx);
plot(orderx, fitted_y, 'r');

% Calculating and displaying the slope of the line
slope = coefficients(1);
text(orderx(end), fitted_y(end), sprintf('Slope: %.2f', slope), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

xlim([1.5,10.5]);
ytickformat('percentage')

xlim([1.5,10.5]);

xlabel('Order')
ylabel('Proportion of weights greater than 1 in each order')
title('Whether higher orders are more likely to collaborate across communities')

% Adding legend
legend('','Original model: Increased probability of cross-community collaboration at higher orders','','Random model: Completely random cross-community probability','','APS 1999-2009');

hold off;
