function [data,n_vec,m_vec,t_vec] = convergence_rate_analysis(varargin)
% --- convergence_rate_analysis() -----------------------------------------
% Numerical evaluation of the convergence rate of GEVO-CI.
%
% 2023-10-30 Robin Forsling

set_latex_interpreter;


% --- USE SAVED DATA ---
if nargin > 0
    v = varargin{1};
    if isstring(v) || ischar(v); load(v);
    else; error('input must be a string or char...')
    end
    M = size(data{1,1},2);


% --- SIMULATE ---
else
    M = 10000;
    n_vec = [6 9];
    m_vec = [1 2 3 4];                                                      
    t_vec = [0.001 0.0001];   
    data = cell(length(t_vec),length(n_vec)); 
    for f = 1:length(t_vec)
        J_thres = t_vec(f);
        for g = 1:length(n_vec)
            n = n_vec(g);
            num_iter = zeros(length(m_vec),M);
            for i = 1:length(m_vec)
                m = m_vec(i);
                for j = 1:M
                    R1 = get_random_covariance(n);
                    R2 = get_random_covariance(n);
                    while min(eig(R2-R1)) >= 0; R1 = get_rnd_cov(n); end                
                    num_iter(i,j) = gevoci_dr(R1,R2,m,J_thres);
                end
            end
            data{f,g} = num_iter;
        end
    end
end


% --- EVALUATION ----------------------------------------------------------
c = get_thesis_colors; clr = cell(6,1); lw = 1.5;
clr{1} = c.red; clr{2} = c.yellow; clr{3} = c.green; 
clr{4} = c.blue; clr{5} = c.purple; clr{6} = c.cyan;
nrow = 1; ncol = length(n_vec)*length(t_vec);

figure(2); clf; 
plot_cnt = 1;
for f = 1:length(t_vec)
    J_thres = t_vec(f);
    for g = 1:length(n_vec)
        subplot(nrow,ncol,plot_cnt); 
        n = n_vec(g); 
        num_iter = data{f,g};

        for i = 1:length(m_vec)
            [v,h] = integer_histogram(num_iter(i,:)); 
            hh = semilogy(v,h/M,'-','DisplayName',sprintf('$m=%d$',m_vec(i))); hh.Color = clr{i}; hh.LineWidth = lw;
            hold on
        end

        title(sprintf('$n=%d$, $\\epsilon=%1.4f$',n,J_thres),'interpreter','latex')
        if true; xlabel('Number of iterations'); end % f == nrow
        if plot_cnt == 1; ylabel('Relative frequency'); end
        
        xlim([0 20])
        legend show; grid on
        plot_cnt = plot_cnt + 1;
    end
end

fprintf('\n\n')
for f = 1:length(t_vec)
    J_thres = t_vec(f);
    for g = 1:length(n_vec)
        n = n_vec(g); 
        num_iter = data{f,g};
        fprintf('J_thres = %1.5f, n = %d\n',J_thres,n)
        for i = 1:length(m_vec)
            m = m_vec(i);
            [typ,ave,stdev] = convergence_rate_stats(num_iter(i,:));  
            fprintf('%d & %d & %2.3f & %2.3f \\\\ \n',m,typ,ave,stdev)
        end
        fprintf('\n')
    end
end

set_fontsize_all(14)


end



