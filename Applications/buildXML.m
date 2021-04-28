function [seq_return] = buildXML(rng_val, name, file_directory, nr_samples, sampleTimes,...
    method, ignoreEnds, fixrate, min_cov, nr_reps, exclude, from, to, recombinationBreakpoints, include)

rng(rng_val);
cov = 1-min_cov;

cov_values = zeros(0,0);

if method > 0 && method < 3
    if method==1
        files = dir([file_directory '/results/*.fasta']);
        d = fastaread([files(3).folder '/' files(3).name]);
    elseif method==2
        files = dir([file_directory '/results/*full.fasta']);
        d = fastaread([files(3).folder '/' files(1).name]);
    end

    c = 1;
    for i = 1 : length(d)
        if to==-1
            d_tmp(i).Sequence = d(i).Sequence(from:end-100);
        else
            d_tmp(i).Sequence = d(i).Sequence(from:to);
        end

        d_tmp(i).Header = d(i).Header;
        
        if getCov(d_tmp(i).Sequence, d_tmp(i).Header)<cov || method==1
            cov_values(c) = getCov(d_tmp(i).Sequence);
            f(c) = d_tmp(i);
            c=c+1;
        end
    end

    if method==1
        g = fopen([file_directory '/data/metadata.tsv']);line = strsplit(fgets(g), '\t');
    elseif method>=2
        tmp = strsplit(files(1).name,'_');
        g = fopen([file_directory '/results/metadata_' tmp{2} '_full.tsv']);line = strsplit(fgets(g), '\t');
    end


    date_id = find(ismember(line, 'date')); c = 1;
    id = cell(0,0);date=cell(0,0);
    while ~feof(g)
        line = strsplit(fgets(g), '\t');
        id{c} = line{1}; date{c} = line{date_id}; c=c+1;
    end
      
    c = 1;    
    for i = 1 : length(f)
        ind = find(ismember(id, f(i).Header));
        tmp2 = strsplit(date{ind}, '-');
        
        if method==1  
            if ~strcmp(date{ind}, '?') && ~ismember(f(i).Header, exclude)
                dat(c) = f(i);
                dat(c).Header = [dat(c).Header '|' tmp2{1}];
                c=c+1;
            end
        elseif method>=2
            if ~strcmp(date{ind}, '?') && ~strcmp(date{ind}, 'NA') && ~contains(f(i).Header, 'ref')
                dat(c) = f(i);
                dat(c).Header = [dat(c).Header '|' date{ind}];
                c=c+1;
            end
        end
    end
    
elseif method==3
    g = fopen([file_directory '/data/example_metadata.tsv']);line = strsplit(fgets(g), '\t');

    date_id = find(ismember(line, 'date')); 
    region_id = find(ismember(line, 'region')); c = 1;
    id = cell(0,0);date=cell(0,0);region=cell(0,0);
    while ~feof(g)
        line = strsplit(fgets(g), '\t');
        id{c} = line{1}; date{c} = line{date_id}; 
        region{c} = line{region_id};c=c+1;
    end
    use_subset = zeros(0,0);    
    uni_reg = unique(region);   
    for i = 1 : length(uni_reg)
        ind = find(ismember(region, uni_reg{i}));
        use_subset = [use_subset, randsample(ind, round(nr_samples/length(uni_reg)))];
    end
    
    for i = 1 : length(include)
         use_subset = [use_subset, find(ismember(id, include{i}))];        
    end
    
    use_subset = unique(use_subset);
    
    id = id(use_subset);
    date = date(use_subset);
    
    
    files = dir([file_directory '/results/split_alignments/*.fasta']);c=1;
    for i = 1: length(files)
        d = fastaread([files(i).folder '/' files(i).name]);
        for j = 1 : length(d)
            if ismember(d(j).Header, id)
                date_str = find(ismember(id, d(j).Header));
                dat(c).Header = [d(j).Header '|' date{date_str}];
                dat(c).Sequence = d(j).Sequence(100:end-100);
                c = c+1;
                fprintf('%d ', c);

            end
        end
    end
    nr_samples = nr_samples*10;
else    
    if method==0
        files = dir([file_directory '/data/*.fasta']);
        d = fastaread([files(1).folder '/' files(1).name]);
        
        
        c = 1;
        for i = 1 : length(d)
            if to==-1
                d_tmp(i).Sequence = d(i).Sequence(from:end-100);
            else
                d_tmp(i).Sequence = d(i).Sequence(from:to);
            end

            d_tmp(i).Header = d(i).Header;
            if getCov(d_tmp(i).Sequence, d_tmp(i).Header)<cov
                cov_values(c) = getCov(d_tmp(i).Sequence);
                f(c) = d_tmp(i);
                c=c+1;
            end
        end

        
        
        c = 1;    
        for i = 1: length(f)
            tmp = strsplit(f(i).Header, '|');
            tmp2 = strsplit(tmp{end}, '-');
            if length(tmp2)==3
                dat(c) = f(i);
                c=c+1;
            end
        end
    elseif method==-1        
        %% read in metadat
        g = fopen([file_directory '/data/samplingrange.csv']);fgets(g);c=1
        id = cell(0,0);date_min=cell(0,0);date_max=cell(0,0);
        while ~feof(g)
            line = strsplit(fgets(g), ',');
            id{c} = line{1}; date_min{c} = line{2};date_max{c} = line{3}; c=c+1;
        end
        fclose(g);       
        
        files = dir([file_directory '/data/*.fasta']);
        f = fastaread([files(1).folder '/' files(1).name]);
        c = 1;    
        for i = 1: length(f)
            tmp = strsplit(f(i).Header, ' ');
            name = strrep(tmp{4}, ',','');
            ind = find(ismember(id, name));
            if ~isempty(ind)
                dat(c) = f(i);
                dat(c).Header = name;
                c=c+1;
            end
        end
    elseif method==-2        
        g = fopen([file_directory '/../data/example_metadata.tsv']);line = strsplit(fgets(g), '\t');
        

        date_id = find(ismember(line, 'date')); c = 1;
        id = cell(0,0);date=cell(0,0);
        while ~feof(g)
            line = strsplit(fgets(g), '\t');
            id{c} = line{1}; date{c} = line{date_id}; c=c+1;
        end

        f = fastaread([file_directory '/global/sample-global.fasta']);

        c = 1;    
        for i = 1: length(f)
            seq = f(i).Header;
            ind = find(ismember(id, seq));
            if ~isempty(ind)
                dat(c) = f(i);
                dat(c).Header = [seq '|' date{ind}]; 
                c=c+1;
            end
        end
    end
end

indices = randsample(1:length(dat), min(length(dat),nr_samples));
data  = dat(indices);
if length(cov_values)>0
    cov_values = cov_values(indices);
end
for reps=0:nr_reps-1

    f = fopen('template/template.xml');
    g = fopen([ 'xmls/' name '_rep' num2str(reps) '.xml'], 'w');

    sample_time = cell(0,0);

    while ~feof(f)
        line = fgets(f);
        if contains(line, 'insert_sequence')
            for i = 1 : length(data)
                if ignoreEnds>0
                    fprintf(g, '\t\t\t\t<sequence id="seq_%s" taxon="%s" totalcount="4" value="%s"/>\n',...
                        data(i).Header, data(i).Header, data(i).Sequence(ignoreEnds:end-ignoreEnds));
                else
                    fprintf(g, '\t\t\t\t<sequence id="seq_%s" taxon="%s" totalcount="4" value="%s"/>\n',...
                        data(i).Header, data(i).Header, data(i).Sequence);
                end
            end

            if sampleTimes
                if ignoreEnds>0
                    fprintf(g, '\t\t\t\t<sequence id="seq_dummy" taxon="dummy" totalcount="4" value="%s"/>\n',...
                        repmat('N',1,length(data(i).Sequence(ignoreEnds:end-ignoreEnds))) );
                else
                    fprintf(g, '\t\t\t\t<sequence id="seq_dummy" taxon="dummy" totalcount="4" value="%s"/>\n',...
                        repmat('N',1,length(data(i).Sequence)) );
                end
            end
        elseif method==1 && contains(line, 'spec="CoalescentWithRecombination"')
            fprintf(g, strrep(line, 'Recombination" maxHeightRatio="1.5"', 'Recombination"  maxHeightRatio="1.1" redFactor="0.1"'));
        elseif contains(line, 'insert_taxa')
            for i = 1 : length(data)
                fprintf(g, '\t\t\t\t<taxon id="%s" spec="Taxon"/>\n', ...
                    data(i).Header);
            end

            if sampleTimes
                fprintf(g, '\t\t\t\t<taxon id="dummy" spec="Taxon"/>\n');

            end

        elseif contains(line, 'insert_sampling_time')
            mrsi = 0;
            for i = 1 : length(data)
                tmp = strsplit(data(i).Header, '|');
                tmp2 = strsplit(tmp{end}, '-');
                mrsi = max(mrsi, str2double(tmp2{1}));

                if method==-1
                    ind = find(ismember(id, data(i).Header));    
                    m = (datenum(date_max{ind}, 'dd-mm-yyyy') - datenum(date_min{ind}, 'dd-mm-yyyy'))/2;
                    sample_time{end+1,1} = data(i).Header;
                    tmp{end} = datestr(m+datenum(date_min{ind}, 'dd-mm-yyyy'), 'yyyy-mm-dd');                
                    mrsi=2003;                
                elseif contains(tmp{end}, 'XX')
                    tmp{end} = strrep(tmp{end}, 'XX-XX','06-01');
                    tmp{end} = strrep(tmp{end}, '-XX','-15');
                    sample_time{end+1,1} = data(i).Header;
                end
                % contains 0 in dates
                if method==1
                    str = tmp{end};
                else
                    str = datestr(datenum(tmp{end}, 'yyyy-mm-dd'),'yyyy-mm-dd');
                end
                    
                if i==length(data) && ~sampleTimes
                    fprintf(g, '\t\t\t\t%s=%s\n', ...
                        data(i).Header, str);
                else
                    fprintf(g, '\t\t\t\t%s=%s,\n', ...
                        data(i).Header, str);
                end
            end

            if sampleTimes
                mrsi = mrsi+1;
                fprintf(g, '\t\t\t\tdummy=%d-01-01\n', mrsi);
            else
                sample_time = cell(0,0);
            end
        elseif contains(line, 'insertBreakPoints')
            fprintf(g, strrep(line, 'insertBreakPoints', sprintf('%d ', recombinationBreakpoints)));
        elseif contains(line, 'insert_prior')       
            fprintf(g, '\t\t\t\t\t\t\t\t\t\t<distribution id="tipprior" spec="recombination.distribution.TipPrior" dateOffset="%d" network="@network">\n', mrsi);
            for i = 1 : length(sample_time)
                fprintf(g, '\t\t\t\t\t\t\t\t\t\t\t<taxonset id="tip.%s" spec="TaxonSet">\n', sample_time{i});
                fprintf(g, '\t\t\t\t\t\t\t\t\t\t\t\t<taxon idref="%s"/>\n', sample_time{i});
                fprintf(g, '\t\t\t\t\t\t\t\t\t\t\t</taxonset>\n');

                if method==-1
                    ind = find(ismember(id, data(i).Header));    
                    min_age = (datenum(date_min{ind}, 'dd-mm-yyyy') - datenum('2003','yyyy'))/365 + 2003;
                    max_age = (datenum(date_max{ind}, 'dd-mm-yyyy') - datenum('2003','yyyy'))/365 + 2003;


                    fprintf(g, '\t\t\t\t\t\t\t\t\t\t\t<distr spec="beast.math.distributions.Uniform" id="Unform.%s" lower="%f" upper="%f"/>\n', sample_time{i}, min_age, max_age);
                else
                    tmp = strsplit(sample_time{i}, '|');
                    tmp2 = strsplit(tmp{end}, '-');

                    fprintf(g, '\t\t\t\t\t\t\t\t\t\t\t<distr spec="beast.math.distributions.Uniform" id="Unform.%s" lower="%s" upper="%d"/>\n', sample_time{i}, tmp2{1}, str2double(tmp2{1})+1);
                end
            end
            fprintf(g, '\t\t\t\t\t\t\t\t\t\t</distribution>\n');
        elseif method==1 && contains(line, 'yyyy-M-dd')
            fprintf(g, strrep(line, 'yyyy-M-dd','yyyy'));
        elseif contains(line, 'insert_operator')
            for i = 1 : length(sample_time)
                fprintf(g, '\t\t\t\t<operator spec="TipReheight" network="@network" size="0.1" weight="0.025">\n');
                fprintf(g, '\t\t\t\t\t<taxonset idref="tip.%s" spec="TaxonSet"/>\n', sample_time{i});
                fprintf(g, '\t\t\t\t</operator>\n');
            end
        elseif contains(line, 'insert_loggers')
            fprintf(g, '\t\t\t\t\t\t<log idref="tipprior"/>\n'); 
%             loci = 0:10:4000
%             for i = 1:length(loci)
%                  fprintf(g, '\t\t\t\t\t\t<log spec="LocusStatsLogger" locus="%d" recombinationNetwork="@network"/>\n', loci(i));
%             end
        elseif fixrate && contains(line, 'mut') 
            if ~contains(line, 'operator') && ~contains(line, 'downParameter')
                fprintf(g, line);
            end
            % do nothing
        else
            fprintf(g, line);
        end
    end
end
if length(cov_values)>0
    ind = find(cov_values==min(cov_values));
    seq_return = data(ind(1));
    seq_return.Sequence = strrep(seq_return.Sequence, '-', 'N');
else
    seq_return = dat(1);
    seq_return.Sequence = strrep(seq_return.Sequence, '-', 'N');
end
end

