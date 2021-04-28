% makes simulation xmls for the coalescent with reassortment
clear

rng(1);

% define the number of repetitions
nr_reps = 10;

% rebuild the xml dirs
system('rm -r inference simulation');
system('mkdir inference simulation');

% define the number of samples
nr_samples = 50;

% define the sampling interval
sampling_interval = 30;

% file that keeps track of the ne and reassortment rates
h = fopen('rates.csv', 'w');fprintf(h, 'run,lambda,Ne,recombination,relRate1,relRate2,relRate3\n');
h2 = fopen('rates_sequences.csv', 'w');fprintf(h2, 'run,lambda,k,m,g,f1,f2,f3,f4\n');

% define params of the lognormal distribution of the Ne
mean_ne = 20;
sigma_ne = 0.25;
mu_ne = log(mean_ne) - sigma_ne^2/2;


% define params of the lognormal distribution of the reassortment rate
mean_mu = 0.0002;
sigma_mu = 0.1;
mu_mu = log(mean_mu) - sigma_mu^2/2;

mean_k = 1;
sigma_k = 1.5;
mu_k = log(mean_k) - sigma_k^2/2;

lambda = [1 2];
for l = 1: length(lambda)
    % make nr reps number of xmls
    for i = 1 : nr_reps
        f_sim = fopen('sim_template.xml');
        
        mean_rea = 0.125/30000/lambda(l);
        sigma_rea = 0.5;
        mu_rea = log(mean_rea) - sigma_rea^2/2;



        % sample the Ne and reassortment rates
        Ne = lognrnd(mu_ne,sigma_ne);
        reassortment = lognrnd(mu_rea, sigma_rea);
        relRate = [0 0 0];


        kappa = lognrnd(mu_k,sigma_k);
        mu = lognrnd(mu_mu, sigma_mu);
        freq = rand(1,4);
        gamma = exprnd(1);


        fprintf(h, '%d,%d,%.12f,%.12f,%.12f,%.12f,%.12f\n', i,lambda(l), Ne, reassortment,relRate(1),relRate(2),relRate(3));
        for j = 1 : length(kappa)
            freq(j,:) = round(freq(j,:)/sum(freq(j,:)),4);
            % ensure that they sum up to 1
            rest = 1;        
            for k = 1 : length(freq(j,:))
                rest = rest-freq(j,k);
            end
            freq(j,end) = freq(j,end)+rest;
            fprintf(h2, '%d,%d,%.12f,%.12f,%.12f%s\n', i,lambda(l), kappa(j), mu(j),gamma(j), sprintf(',%.4f', freq(j,:)));
        end

        % open the simulation xml
        g = fopen(sprintf('simulation/sim_%d.%d.xml', i,lambda(l)), 'w');

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
            elseif ~isempty(strfind(line, 'insert_lambda'))
                fprintf(g, strrep(line, 'insert_lambda', num2str(lambda(l))));
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
        for r = 1 : 3
            % build the inference xml
            f_inf = fopen('inf_template.xml');


            % open the inference xml
            g = fopen(sprintf('inference/inf_hky_%d_rep%d.%d.xml', i, r, lambda(l)), 'w');

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
                    line = strrep(line, 'insert_sim_file_name', sprintf('sim_%d.%d', i, lambda(l)) );
                    fprintf(g, line);
                    segmentcount = segmentcount + 1;
                else
                    fprintf(g, '%s', line);
                end
            end
            fclose(f_inf); fclose(g);
        end
    end
end
fclose(h);
