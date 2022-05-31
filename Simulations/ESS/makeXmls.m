% makes simulation xmls for the coalescent with reassortment
clear

rng(1);

% define the number of repetitions
nr_reps = 300;

% rebuild the xml dirs
system('rm -r inference simulation');
system('mkdir inference simulation');

% define the number of samples
diff_samples = [25 50 100 200];

% define the sampling interval
sampling_interval = 30;

% file that keeps track of the ne and reassortment rates
h = fopen('rates.csv', 'w');fprintf(h, 'run,nr_samples,Ne,recombination,relRate1,relRate2\n');
h2 = fopen('rates_sequences.csv', 'w');fprintf(h2, 'run,m\n');

% define params of the lognormal distribution of the Ne
mean_ne = 10;
sigma_ne = 0.5;
mu_ne = log(mean_ne) - sigma_ne^2/2;

% define params of the lognormal distribution of the reassortment rate
mean_rea = 0.5/30000;
sigma_rea = 0.5;
mu_rea = log(mean_rea) - sigma_rea^2/2;

mean_mu = 0.0008;
sigma_mu = 0.5;
mu_mu = log(mean_mu) - sigma_mu^2/2;

mean_k = 1;
sigma_k = 1.5;
mu_k = log(mean_k) - sigma_k^2/2;

curr_run = 0;
% make nr reps number of xmls
for i = 1 : nr_reps
    nr_samples = randi(190)+10;

    curr_run=curr_run+1;
    f_sim = fopen('sim_template.xml');

    % sample the Ne and reassortment rates
    Ne = lognrnd(mu_ne,sigma_ne);
    reassortment = rand(1)*0.00001+0.00001;
    relRate = [0 0];


    kappa = lognrnd(mu_k,sigma_k);
    mu = lognrnd(mu_mu, sigma_mu);
%         freq = rand(1,4);
%         gamma = exprnd(1);


    fprintf(h, '%d,%d,%.12f,%.12f,%.12f,%.12f\n', curr_run, nr_samples, Ne, reassortment,relRate(1),relRate(2));
    fprintf(h2, '%d,%.12f\n', i, mu);

    % open the simulation xml
    g = fopen(sprintf('simulation/sim_%d.xml', curr_run), 'w');
    clear time

    while ~feof(f_sim)
        line = fgets(f_sim);
        if ~isempty(strfind(line, 'insert_taxa'))
            for j = 1 : nr_samples
                fprintf(g,'\t\t\t\t<taxon spec="Taxon" id="t%d"/>\n', j);
            end             
        elseif ~isempty(strfind(line, 'insert_sampling_times'))
            time = zeros(nr_samples,1);
            for j = 1 : nr_samples
                time(j) = rand*sampling_interval;
                if j==nr_samples
                    fprintf(g,'\t\t\t\tt%d=%f\n', j, time(j));
                else
                    fprintf(g,'\t\t\t\tt%d=%f,\n', j, time(j));
                end
            end 
        elseif ~isempty(strfind(line, 'insert_Ne'))
            fprintf(g, strrep(line, 'insert_Ne', num2str(Ne)));
        elseif ~isempty(strfind(line, 'insert_reassortment'))
            fprintf(g, strrep(line, 'insert_reassortment', num2str(reassortment)));
        elseif ~isempty(strfind(line, 'insert_rel_rates'))
            fprintf(g, strrep(line, 'insert_rel_rates', sprintf('%.12f ', relRate)));
        elseif ~isempty(strfind(line, 'insert_mut'))
            fprintf(g, strrep(line, 'insert_mut', num2str(mu)));
        elseif ~isempty(strfind(line, 'insert_kappa'))
            fprintf(g, strrep(line, 'insert_kappa', num2str(kappa)));
        elseif ~isempty(strfind(line, 'insert_freq'))
            fprintf(g, strrep(line, 'insert_freq', num2str(freq)));
        elseif ~isempty(strfind(line, 'insert_gamma'))
            fprintf(g, strrep(line, 'insert_gamma', num2str(gamma)));

        else
            fprintf(g, '%s', line);
        end
    end
    fclose(f_sim); fclose(g);

    % make 3 replicates
    for r = 1 : 1
        % build the inference xml
        f_inf = fopen('inf_template.xml');


        % open the inference xml
        g = fopen(sprintf('inference/inf_hky_%d_rep%d.xml', curr_run, r), 'w');

        % keep track of the segment count for the nexus file name
        segmentcount = 1;segmentcount2=1;

        while ~feof(f_inf)
            line = fgets(f_inf);
            if ~isempty(strfind(line, 'insert_taxa'))
                for j = 1 : nr_samples
                    fprintf(g,'\t\t\t\t<taxon spec="Taxon" id="t%d"/>\n', j);
                end             
            elseif ~isempty(strfind(line, 'insert_sampling_times'))
                for j = 1 : nr_samples
                    if j==nr_samples
                        fprintf(g,'\t\t\t\tt%d=%f\n', j, time(j));
                    else
                        fprintf(g,'\t\t\t\tt%d=%f,\n', j, time(j));
                    end
                end 
            elseif ~isempty(strfind(line, 'initial_Ne'))
                fprintf(g, strrep(line, 'initial_Ne', num2str(exprnd(1))));
            elseif ~isempty(strfind(line, 'initial_rel_rates'))
                fprintf(g, strrep(line, 'initial_rel_rates', sprintf('%.3f ',normrnd(0,1,1,2) )));
            elseif ~isempty(strfind(line, 'initial_reassortment'))
                fprintf(g, strrep(line, 'initial_reassortment', num2str(exprnd(mean_rea))));
            elseif ~isempty(strfind(line, 'insert_sim_file_name'))
                line = strrep(line, 'insert_sim_file_name', sprintf('sim_%d', curr_run) );
                fprintf(g, line);
                segmentcount = segmentcount + 1;
            else
                fprintf(g, '%s', line);
            end
        end
        fclose(f_inf); fclose(g);
    end
end

fclose(h);