%% Computation

clearvars -except k;
cd("/Users/liuyuyan/Documents/MATLAB/projects/sci_coop_nw/");
load("/Users/liuyuyan/Documents/MATLAB/projects/sci_coop_nw/data/nw_a/bip_year.mat");
bip = readmatrix("./data/nw_a/authorid_paperid.csv");
bip = bip + 1; % Adjust IDs to MATLAB's 1-based indexing
tem = bip(:,2);
bip(:,2) = bip(:,1);
bip(:,1) = tem; % Swap columns: paper-author -> author-paper
clear tem;

year_idx = find(year2 >= 1999); % Control year
bip = bip(year_idx,:);

paper_frequence = tabulate(bip(:,1));
order_frequence = tabulate(paper_frequence(:,2));
sorted_order_frequence = sortrows(order_frequence, 1);
sorted_order_frequence(1,:) = [];

for i = 1:length(sorted_order_frequence(:,1))
    order = sorted_order_frequence(i,1);
    sorted_order_edge_frequence(i,1) = sorted_order_frequence(i,2) * order * (order - 1) / 2;
end
sorted_order_edge_frequence = sorted_order_edge_frequence .* 100 ./ sum(sorted_order_edge_frequence);

%% Plotting
plot(sorted_order_frequence(:,1), sorted_order_edge_frequence(:) ./ sorted_order_frequence(:,3));

%% Color Sets for Plotting
color1 = [142, 207, 201] / 255; % #8ECFC9
color2 = [255, 190, 122] / 255; % #FFBE7A
color3 = [250, 127, 111] / 255; % #FA7F6F
color4 = [130, 176, 210] / 255; % #82B0D2
color5 = [190, 184, 220] / 255; % #BEB8DC
color6 = [231, 218, 210] / 255; % #E7DAD2

mycolorset6 = [color1; color2; color3; color4; color5; color6];
mycolorset3 = [77 133 189; 247 144 61; 89 169 90] / 255;
mycolorset2 = [[238/255, 136/255, 129/255]; [141/255, 183/255, 219/255]];
mycolorset4 = [243 175 171; 243 204 91; 158 193 222; 147 205 209] / 255;

%% Calculating the Largest Connected Group of Authors and Papers
[all_in_one] = bip2one(bip); % Weighted all_in_one
[C, sizes] = conncomp(graph(all_in_one));
[~, idx] = max(sizes);
L = all_in_one(C == idx, C == idx);

%% Bar Graph - Data for the First 20 Orders with Color Modification
figure
height = gcf().Position(4); % Get current figure height
width = height * 2.4; % Calculate required width based on height
set(gcf(), 'Position', [gcf().Position(1), gcf().Position(2), width, height]); % Set new figure dimensions

bar(sorted_order_frequence(:,1), [sorted_order_frequence(:,3), sorted_order_edge_frequence(:)], 'grouped', 'BarWidth', 1.7);
ytickformat('percentage'); % Format y-axis as percentage

% Set bar colors
h = gca;
h.Children(1).FaceColor = mycolorset2(1,:); % Color of the first bar
h.Children(2).FaceColor = mycolorset2(2,:); % Color of the second bar

% Add chart title and axis labels
xlabel('Order');
ylabel('Percentage');
xlim([0, 20]);
legend('Number of Papers', 'Number of Edges');
set(gca(), 'FontSize', 23);
set(gca(), 'LineWidth', 2);

%% Bar Graph for Orders Higher than 20
figure;
height = gcf().Position(4); % Get current figure height
width = height * 3.5; % Calculate required width based on height
set(gcf(), 'Position', [gcf().Position(1), gcf().Position(2), width, height]); % Set new figure dimensions

% Right y-axis for edges data
yyaxis right;
b2 = bar(sorted_order_frequence(20:end,1), sorted_order_edge_frequence(20:end), 'FaceColor', mycolorset2(1,:), 'EdgeColor', 'none', 'BarWidth', 3);
ylabel('Percentage of Edges');
ytickformat('percentage');
ylim([0 10]);
set(gca, 'YColor', [0 0 0]);

% Left y-axis for papers data
yyaxis left;
b1 = bar(sorted_order_frequence(20:end,1), sorted_order_frequence(20:end,3), 'FaceColor', mycolorset2(2,:), 'EdgeColor', 'none', 'BarWidth', 3);
ylabel('Percentage of Papers');
ytickformat('percentage');
ylim([0 0.01]);
title('Higher than 20 Order');
set(gca, 'YColor', [0 0 0]);
set(gca(), 'FontSize', 35);
set(gca(), 'LineWidth', 2.8);
xticks([20, 100, 200, 300, 400, 500, 600]);
xticklabels([20, 100, 200, 300, 400, 500, 600]);
print('2.1-a2.2.eps', '-depsc');
