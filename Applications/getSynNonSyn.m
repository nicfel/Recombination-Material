% get synonynmous and non-syn differences on different parts of the genome
% of sars-like viruses
clear


name = {'mers',...
    'sars-like',...
    '229e',...
    'nl63',...
    'oc43'};

file_directory = {'mers',...
    '../../sars-like-cov',...
    '../../seasonal-cov/229e',...
    '../../seasonal-cov/nl63',...
    '../../seasonal-cov/oc43'};

method = [0, 1, 2, 2, 2];




for i = 4 : length(name)-1
    disp(name{i})
    if method(i)==0
        files = dir([file_directory{i} '/data/*.fasta']);
        d = fastaread([files(1).folder '/' files(1).name]);
    elseif method(i)==1
        files = dir([file_directory{i} '/results/*.fasta']);
        d = fastaread([files(3).folder '/' files(3).name]);
    elseif method(i)==2
        files = dir([file_directory{i} '/results/*full.fasta']);
        d = fastaread([files(3).folder '/' files(1).name]);
    end
    
    protein = tdfread(['reference/' name{i} '.tsv'],'\t');
    gb = genbankread(['reference/' name{i} '.gb']);

    f = fopen(['config/' name{i} '.tsv'], 'w');
    fprintf(f, 'protein\tsynonymous\tmutations\tlength\n');
        
    length_values = zeros(0,3);
        
    for j = 1 : length(gb.CDS)
        disp(j)
        % find stat position
        seq = upper(gb.Sequence(gb.CDS(j).indices(1):gb.CDS(j).indices(1)+15));
        for k = 1 : length(d)
            start = strfind(d(k).Sequence, seq);            
            if ~isempty(start)
                break;
            end
        end    
        clear prot seq 
        c=1;
        for k = 1 : length(d)
            seq_tmp = d(k).Sequence(start:start+(gb.CDS(j).indices(end)-gb.CDS(j).indices(1)));
            seq_tmp(seq_tmp=='-')='N';
            if sum(seq_tmp=='N')/length(seq_tmp) < 0.1
                seq(c,:) = seq_tmp;
                prot(c,:) = nt2aa(seq_tmp, 'ACGTOnly', false);
                c=c+1;
            end                
        end
        
        count_seq = 0;
        for k = 1 : size(seq,2)
            u = unique(seq(:,k));
            u(u=='N') = [];
            count_seq = count_seq + length(u)-1;
        end
        
        count_prot = 0;
        for k = 1 : size(prot,2)
            u = unique(prot(:,k));
            u(u=='X') = [];
            count_prot = count_prot + length(u)-1;
        end
        
        length_values(j,:) = [count_prot, count_seq, gb.CDS(j).indices(end)-gb.CDS(j).indices(1)];
        fprintf(f, '%s\t%d\t%d\t%d\n', gb.CDS(j).gene, count_prot, count_seq, gb.CDS(j).indices(end)-gb.CDS(j).indices(1));       
    end   
    fclose(f);
    
    if i==1 || i==3
    
        length_ratios(1) = sum(length_values(1,1))/(sum(length_values(1,2))-sum(length_values(1,1)));
        length_ratios(2) = sum(length_values(3,1))/(sum(length_values(3,2))-sum(length_values(3,1)));
        length_ratios(3) = sum(length_values(4:end,1))/(sum(length_values(4:end,2))-sum(length_values(4:end,1)));
    elseif i==2
        length_ratios(1) = sum(length_values(1:2,1))/(sum(length_values(1:2,2))-sum(length_values(1:2,1)));
        length_ratios(2) = sum(length_values(3,1))/(sum(length_values(3,2))-sum(length_values(3,1)));
        length_ratios(3) = sum(length_values(4:end,1))/(sum(length_values(4:end,2))-sum(length_values(4:end,1)));
    elseif i==4
        length_ratios(1) = sum(length_values(1,1))/(sum(length_values(1,2))-sum(length_values(1,1)));
        length_ratios(2) = sum(length_values(2,1))/(sum(length_values(2,2))-sum(length_values(2,1)));
        length_ratios(3) = sum(length_values(3:end,1))/(sum(length_values(3:end,2))-sum(length_values(3:end,1)));
    end
    plot(length_values(:,1)./length_values(:,3)); hold on
    plot((length_values(:,2) -length_values(:,1))./length_values(:,3))


    avg_ratio = sum(length_values(1:end,1))/(sum(length_values(1:end,2))-sum(length_values(1:end,1)));
%     scatter(1:3,length_ratios/avg_ratio); hold on
%     title(name{i})

end
