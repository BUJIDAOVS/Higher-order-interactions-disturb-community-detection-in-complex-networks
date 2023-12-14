sim_big_author = [];
for n = 1:50

    little_auth = n;

    %% Community ID-PACS codes for networks less than 15 (excluding minor authors)
    com_num_less_then_15 = size(unique(S_less_then_15),1);
    com_pacs_less_then_15 = cell(com_num_less_then_15, 1);
    for i = 1:length(com_pacs_less_then_15)
        com_pacs_less_then_15{i} = [];
    end

    % Assigning community IDs to authors
    for i = [1:size(S_less_then_15,1)]
        if length(find(apyp(:,1)==i)) >= little_auth
            tem = apyp(find(apyp(:,1)==i),4:end);
            com_pacs_less_then_15{S_less_then_15(i)} = [com_pacs_less_then_15{S_less_then_15(i)},tem(:)'];
            com_pacs_less_then_15{S_less_then_15(i)}(com_pacs_less_then_15{S_less_then_15(i)}==0)=[];
            clear tem;
        end
    end

    %% Community ID-PACS codes for all orders (excluding minor authors)
    com_num = size(unique(S),1);
    com_pacs = cell(com_num, 1);
    for i = 1:length(com_pacs)
        com_pacs{i} = [];
    end

    % Assigning community IDs to authors
    for i = [1:size(S,1)]
        if length(find(apyp(:,1)==i)) >= little_auth
            tem = apyp(find(apyp(:,1)==i),4:end);
            com_pacs{S(i)} = [com_pacs{S(i)},tem(:)'];
            com_pacs{S(i)}(com_pacs{S(i)}==0)=[];
            clear tem;
        end
    end

    %% Community author field similarity (calculated using networks of orders less than or equal to 15)
    clear sim;

    % Creating a cell array
    C = com_pacs_less_then_15;
    % Or C = com_pacs;

    % Using cellfun function to get logical indices of non-empty elements
    idx = cellfun(@isempty, C);
    idx = ~idx;

    % Calculating frequency distribution
    fre = tabulate(S_less_then_15);
    com_id_useful = fre((find(fre(:,2) > 1)),1);
    % Or com_id_useful = find(idx);

    for j = 1:length(com_id_useful)
        clear tem_fre i;
        i = com_id_useful(j);
        if isempty(com_pacs_less_then_15{i})
            continue
        end
        tem_fre = tabulate(fix(com_pacs_less_then_15{i}));
        sorted = sort(tem_fre(:,3), 'descend');
        sim(i) = sum(sorted(1:1))/100;
    end
    sim(sim == 0) = [];
    length_sim_less15 = length(sim);
    mean_sim_less15 = mean(sim);
    %% Community author field similarity (calculated using networks of all orders)
    clear sim;

    % Creating a cell array
    % Or C = com_pacs_less_then_15;
    C = com_pacs;

    % Using cellfun function to get logical indices of non-empty elements
    idx = cellfun(@isempty, C);
    idx = ~idx;

    % Calculating frequency distribution
    fre = tabulate(S);
    com_id_useful = fre((find(fre(:,2) > 1)),1);  % Excluding all single-member communities
    % Or com_id_useful = find(idx);   % Keeping single-member communities with PACS codes
    for j = 1:length(com_id_useful)
        clear tem_fre i;
        i = com_id_useful(j);
        if isempty(com_pacs{i})
            continue
        end
        tem_fre = tabulate(fix(com_pacs{i}));
        sorted = sort(tem_fre(:,3), 'descend');
        sim(i) = sum(sorted(1:1))/100;
    end
    sim(sim == 0) = [];
    length_sim_all = length(sim);
    mean_sim_all = mean(sim)

    sim_big_author(n,:) = [mean_sim_less15,length_sim_less15,mean_sim_all,length_sim_all];
end


%% Data: data/nw_a/sim_big_author.mat
% The calculation part has been completed previously, now just load the data for plotting
load("data/nw_a/sim_big_author.mat")


%% 
figure
color1 = [142, 207, 201] / 255;  % #8ECFC9
color2 = [255, 190, 122] / 255;  % #FFBE7A
color3 = [250, 127, 111] / 255;  % #FA7F6F
color4 = [130, 176, 210] / 255;  % #82B0D2
color5 = [190, 184, 220] / 255;  % #BEB8DC
color6 = [231, 218, 210] / 255;  % #E7DAD2

mycolorset6 = [color1;color2;color3;color4;color5;color6];



mycolorset3 = [77 133 189;247 144 61;89 169 90]/255;

mycolorset2 = [[238/255, 136/255, 129/255];[141/255, 183/255, 219/255]];

mycolorset4 = [243	175	171	;243	204	91	;158	193	222	;147	205	209	]/255


%lineProperties = {'-o', 'Color', mycolorset2(1,:),'MarkerSize', 4, 'markerfacecolor', mycolorset2(1,:), 'MarkerEdgeColor', 'none', 'LineWidth', 9};

plot([1:50],sim_big_author(:,1), '-o', 'Color', mycolorset2(2,:),'MarkerSize', 9, 'markerfacecolor', mycolorset2(2,:), 'MarkerEdgeColor', 'none', 'LineWidth', 3);
hold on
plot([1:50],sim_big_author(:,3), '-s', 'Color', mycolorset2(1,:),'MarkerSize', 9, 'markerfacecolor', mycolorset2(1,:), 'MarkerEdgeColor', 'none', 'LineWidth', 3);
% title('Consistency of field for community authors after excluding minor authors')
xlabel('Publication Threshold')
ylabel('Consistency')

legend('Low-order Network (<15)','Complete Network')
% Increase font size of x and y axis labels, specify value for boldness
set(gca(), 'FontSize', 23)

% Turn off the border of the legend
legend('boxoff','FontSize', 19, 'Location', 'northwest')

xlim([0,51]);


% Set line width of the border
set(gca(), 'LineWidth', 2)
print('3.3-a.eps','-depsc', '-r600')

%% 
figure
plot([1:50],sim_big_author(:,2), '-o', 'Color', mycolorset2(2,:),'MarkerSize', 9, 'markerfacecolor', mycolorset2(2,:), 'MarkerEdgeColor', 'none', 'LineWidth', 3);
hold on
plot([1:50],sim_big_author(:,4), '-s', 'Color', mycolorset2(1,:),'MarkerSize', 9, 'markerfacecolor', mycolorset2(1,:), 'MarkerEdgeColor', 'none', 'LineWidth', 3);
% title('Number of non-empty communities after excluding minor authors')
xlabel('Publication Threshold')
ylabel('Valid Communities')

legend('Low-order Network (<15)','Complete Network')
% Increase font size of x and y axis labels, specify value for boldness
set(gca(), 'FontSize', 23)

% Turn off the border of the legend
legend('boxoff','FontSize', 19, 'Location', 'northwest')

xlim([0,51]);


% Set line width of the border
set(gca(), 'LineWidth', 2)
print('3.3-b.eps','-depsc', '-r600')



