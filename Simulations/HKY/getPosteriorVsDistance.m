% analyses the detection capabilities of the reassortment events
clear
% get all tree files
tree_files = dir('target/*.trees');

% file to print the results to 
g = fopen('event_support.csv', 'w');
fprintf(g, 'dist,support,nrSegments,information\n');


for i = 1:length(tree_files)
    disp(i)
    % read in the network
    f = fopen(['target/' tree_files(i).name]);
    while ~feof(f)
        line = fgets(f);
        if length(line)>1000
            % get the annotated nodes
            meta_data = regexp(line, '\[(.*?)\]', 'match');
        end
    end
    fclose(f);
    
    % read in the reassortment distances 
    f = fopen(['target/' strrep(tree_files(i).name, '.trees', '.txt')]);
    
    clear height min_dist   
    while ~feof(f)
        line = strsplit(strtrim(fgets(f)), '\t');
        for j = 1 : length(line)
            tmp = strsplit(line{j}, ',');
            min_par = 1000000000;
            height(j) = 0;
            for k = 2 : length(tmp)
                tmp2 = strsplit(tmp{k}, ':');
                height(j) = str2double(tmp2{2});
                parent_height = str2double(tmp2{3});
                if parent_height < min_par
                    min_par = parent_height;
                end
            end
            min_dist(j) = min_par - height(j);
        end
        
    end
    fclose(f);
    
    found = false(size(height));
    
    posterior = ones(size(height))*-1;
    
    % match the nodes from the meta_data to distances
    for j = 1 : length(meta_data)
        % check if the event is either a coalescent or reassortment event
        if contains(meta_data{j}, 'targetHeight')
            % check if there is a corresponding node 
            tmp = strsplit(meta_data{j}, '=');
            tmp2 = strsplit(tmp{9}, ',');
            ind = find(abs(height-str2double(tmp2{1}))<0.00000001);
            if ~isempty(ind)
                tmp3 = strsplit(tmp{8}, ',');
                posterior(ind) = str2double(tmp3{1});
            end
        end
    end
    
    % 
    for j = 1 : length(posterior)       
        tmp = strsplit(tree_files(i).name, '.');
        if min_dist(j)<100000
            if contains(tree_files(i).name, 'high')
                information='high';
            else
                information='low';
            end

            if contains(tree_files(i).name, '2seg')
                fprintf(g, '%f,%f,%d,%s\n', min_dist(j), posterior(j), 2, information); 
            elseif contains(tree_files(i).name, '3seg')
                fprintf(g, '%f,%f,%d,%s\n', min_dist(j), posterior(j), 3, information); 
            else
                fprintf(g, '%f,%f,%d,%s\n', min_dist(j), posterior(j), 4, information);
            end
        end
    end
end
fclose(g);