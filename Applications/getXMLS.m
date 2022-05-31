clear

name = {'mers',...
    'sars-like',...
    '229e',...
    'nl63',...
    'oc43',...
    'sars2-like',...
    'sars-cov2'};

file_directory = {'mers',...
    '../../sars-like-cov',...
    '../../seasonal-cov/229e',...
    '../../seasonal-cov/nl63',...
    '../../seasonal-cov/oc43',...
    '../../sars-like-cov',...
    '../../ncov'};

methods = [0, 1, 2, 2, 2, 1, 3];
nr_samples = [100, 100, 100, 100, 100, 100, 100];
min_cov = [-1, -1, 0.1, 0.1, 0.2, -1, 1];
% samplesamples = [false, false, true, true, true, false, false];
samplesamples = [false, false, false, false, false, false, false];

fixrate = [false, true, false, false, false, true, false];

rngs = [0654, 68741, 387421, 3543, 534453, 564555,  4685312];


nr_reps = 3;

exclude = {'BJ02', 'CUHK_AG03','TW6', 'GZ0402','GZ0401', 'HC/GZ/32/03', 'PC4_136','PC4_13',...
    'Rs4231', 'Rs4247','Rs4081','Rs4874','Rs4092','Rs7327','Rs9401','Rs4255','Rs3367','Rs4084',...
    'Rf4092', 'As6526', 'RsSHC014', 'WIV1', 'Sin842', 'civet020', 'Rs672', 'BtKY72', 'bat_SL_CoVZXC21', 'bat_SL_CoVZC45',...
    'LYRa11', 'Wuhan/WIV05/2019', 'USA/CA1/2020', 'Sweden/01/2020','USA/IL1/2020',...
    'Japan/KY-V-029/2020', 'Nepal/61/2020'};

trim_ends = 100;

%% Makes xml for spike protein only

map = tdfread('reference/name_map.tsv','\t');

for i = 1 : length(name)-1
    if i < 6
        disp(name{i})
        protein = tdfread(['reference/' name{i} '.tsv'],'\t');
        gb = genbankread(['reference/' name{i} '.gb']);

        index = -1;
        for j = 1 : size(protein.protein,1)
    %         disp(protein.protein(j,:))
            for k = 1 : size(map.protein,1)
                if contains(map.protein(k,:), protein.protein(j,:))
                    if contains(map.to(k,:), 'Spike')
                        index = j;
                    end
                end
            end
        end

        from = protein.from(index);
        to = protein.to(index);
        buildXML(rngs(i), [name{i} '_all'], file_directory{i}, nr_samples(i), samplesamples(i), methods(i), 0, fixrate(i), 0.1, nr_reps, {''}, 100, -1, [from, to]);    
    else
        buildXML(rngs(i), [name{i} '_all'], file_directory{i}, nr_samples(i), samplesamples(i), methods(i), 0, fixrate(i), 0.9, nr_reps, exclude, 100, -1, [1]);    
    end
end

%%
% i = length(name);
% % VanInsberghe recombinants
% include = {'USA/CA-CZB-1437/2020',...
%     'Beijing/DT-travelIT05/2020',...
%     'USA/CA-CSMC109/2020',...
%     'Wuhan/HB-WH6-246/2020',...
%     'England/201090235/2020'};

% buildXML(rngs(i), [name{i} '_recomb'], file_directory{i}, nr_samples(i), samplesamples(i), methods(i), 0, fixrate(i), 0.25, nr_reps, exclude, 100, -1, include);
% buildXML(rngs(i), [name{i} '_all'], file_directory{i}, nr_samples(i), samplesamples(i), methods(i), 0, fixrate(i), 0.25, nr_reps, exclude, 100, -1, '');
% 
% include = {'England/20118111104/2020',...
%     'SouthAfrica/N00380/2020',...
%     'SouthAfrica/NHLS-UCT-GS-0808/2020',...
%     'England/CAMC-B209F5/2020',...
%     'England/QEUH-BA2A70/2020',};
% buildXML(rngs(i), [name{i} '_N501'], '../../ncov_recomb', 50, samplesamples(i), methods(i), 0, fixrate(i), 0.25, nr_reps, exclude, 100, -1, '');


%%
for i = 1 : 5
    disp(name{i})
    protein = tdfread(['reference/' name{i} '.tsv'],'\t');
    gb = fastaread(['reference/' name{i} '.fasta']);
    
    ind_s1 = zeros(0,0);
    ind_s2 = zeros(0,0);
    
    if i>=3
    
        f = fopen(['../../seasonal-cov/' name{i} '/config/' name{i} '_s1_reference.gb']);

        clear prot_length_s1 prot_length_s2 ind_s1 ind_s2

        % find the starting points for s1 and s2
        while ~feof(f)
            line = strtrim(fgets(f));
            if contains(line, '..')
                prot_length_tmp = strsplit(line, '..');
                prot_length_s1 = str2double(prot_length_tmp{end});
                if ~contains(prot_length_tmp{1}, '1')
                    error('should be a 1')
                end
            elseif contains(line, 'ORIGIN')
                line = fgets(f);
                tmp = strsplit(strtrim(strrep(strtrim(line), '1', '')));
                ind_s1 = strfind(gb(2).Sequence, upper(tmp{1}));
                disp(upper(tmp{1}))
                break;
            end
        end
        fclose(f);

        f = fopen(['../../seasonal-cov/' name{i} '/config/' name{i} '_s2_reference.gb']);
        % find the starting points for s1 and s2
        while ~feof(f)
            line = strtrim(fgets(f));
            if contains(line, '..')
                prot_length_tmp = strsplit(line, '..');
                prot_length_s2 = str2double(prot_length_tmp{end});
                if ~contains(prot_length_tmp{1}, '1')
                    error('should be a 1')
                end
            elseif contains(line, 'ORIGIN')
                line = fgets(f);
                tmp = strsplit(strtrim(strrep(strtrim(line), '1', '')));
                ind_s2 = strfind(gb(2).Sequence, upper(tmp{1}));
                disp(upper(tmp{1}))
                break;
            end
        end
        fclose(f);
    else
        disp(name{i})
        index = -1;
        for j = 1 : size(protein.protein,1)
    %         disp(protein.protein(j,:))
            for k = 1 : size(map.protein,1)
                if contains(map.protein(k,:), protein.protein(j,:))
                    if contains(map.to(k,:), 'Spike')
                        index = j;
                    end
                end
            end
        end

        from = protein.from(index);
        to = protein.to(index);

        
        f = fopen('./config/subunit.tsv');fgets(f);
        
        if i==1
            files = dir([file_directory{i} '/data/*.fasta']);
            d = fastaread([files(1).folder '/' files(1).name]);
        elseif i==2
            files = dir([file_directory{i} '/results/*.fasta']);
            d = fastaread([files(3).folder '/' files(3).name]);
        end

        while ~feof(f)
            line = strsplit(strtrim(fgets(f)));
            if strcmp(line{1}, name{i})
                if strcmp(line{2}, 's1')
                    use_seq_nr=1;
                    while isempty(ind_s1)
                        ind_s1 = strfind(d(use_seq_nr).Sequence, line{3});
                        disp(ind_s1)
                        use_seq_nr = use_seq_nr+1;
                    end
                else
                    use_seq_nr=1;
                    while isempty(ind_s2)
                        ind_s2 = strfind(d(use_seq_nr).Sequence, line{3});
                        use_seq_nr = use_seq_nr+1;
                    end
                    prot_length_s2 = to-ind_s2;
                    prot_length_s1 = ind_s2-ind_s1;
                    
                end
            end
        end
        fclose(f);
    end
    
%     index = -1;
%     for j = 1 : size(protein.protein,1)
% %         disp(protein.protein(j,:))
%         for k = 1 : size(map.protein,1)
%             if contains(map.protein(k,:), protein.protein(j,:))
%                 if contains(map.to(k,:), 'Spike')
%                     index = j;
%                 end
%             end
%         end
%     end
    from = ind_s1;
    to = ind_s2 + prot_length_s2 ;    
    bP = prot_length_s1;
    if i==5
        disp('lla')
    end
    buildXML(rngs(i), [name{i} '_subunits'], file_directory{i}, nr_samples(i), samplesamples(i), methods(i), 0, fixrate(i), 0.25, nr_reps, {''}, from, to, bP);
end

