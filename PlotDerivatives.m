% Partial derivatives
% 
% mikael.mieskolainen@cern.ch, 2019
clear; close all; clc;

addpath src

%% Loop over
for N = [2:5]
    
    all = {};
    
    % Derivative order
    for order = 0:1:2
        
        fprintf('N = %d, Order = %d \n', N, order);
        
        % Symbolic variables
        A = sym(amat(N));
        p = sym('p',[2^N-1 1]);
        u = sym('u');
        
        %% Get partial derivatives
        hat = partialderivative(N);
        
        y = inv(A) *( exp(-u*A*p) - 1) / (exp(-u)-1);
        hat = diff(y, u, order);

        % Maximum entropy input (uniform) + noise for visualization
        pval = ones(1, 2^N-1)+rand(1,2^N-1)*0.1; pval = pval / sum(pval);
        
        % Create symbolic variables
        for i = 1:2^N-1
            eval(sprintf('p%d = sym(''p%d'');', i, i)); 
        end
        
        %
        % mu-values
        uval = linspace(-20.01, 20.01, 200)+1e-2; % to avoid division by zero at u = 0
        
        % Create symbolic substitution expression
        str = 'val = subs(hat, [';
        for i = 1:2^N-1
            str = [str, sprintf('p%d ', i)];
        end
        str = [str, 'u], [pval uval(i)]);'];
        vals = zeros(2^N-1,length(uval));
        for i = 1:length(uval)
            fprintf('N = %d, uval = %d/%d \n', N, i, length(uval));
            eval(str);
            vals(:,i) = double(val);
        end
        
        % Plot
        f1 = figure;
        
        ylimit = max(abs([min(vals(:))*1.1 max(vals(:))*1.1]));
        
        % Rectangle
        rectangle('Position',[-19.9 -ylimit 20 2*ylimit], 'FaceColor', [ones(3,1)*0.96], 'Edgecolor', ones(3,1)); hold on;
        
        % Function
        h1 = plot(linspace(min(uval), max(uval), 3), zeros(3,1), 'k-.'); hold on;
        h2 = plot(uval, vals');
        xlabel('$\mu$','interpreter','latex');
        xticks([-20:5:20]);
        axis([min(uval) max(uval) -ylimit ylimit]);
        
        set(gca, 'Layer', 'top'); grid off; % axis on top
        
        all{order+1} = vals;
        
        legends = {};
        for c = 1:length(y)
            legends{end+1} = sprintf('$y_{%d}$',c);
        end
        title(sprintf('$N = %d$', N),'interpreter','latex');
        
        if (N <= 4)
            l = legend(h2, legends);
            set(l,'interpreter','latex', 'location','southeast');
            
            % with transparency
            transparency = 0.85;
            if (N >= 4)
                set(l.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[ones(3,1); transparency]));
            else
                set(l,'box','off');    
            end
        end
        
        if     (order == 0)
            ylabel('$\mathbf{y}$','interpreter','latex');
        elseif (order == 1)
            ylabel('$\partial \mathbf{y}/ \partial \mu$','interpreter','latex');
        else
            ylabel(sprintf('$\\partial^{%d} \\mathbf{y}/ \\partial \\mu^{%d}$', order, order), 'interpreter', 'latex');
        end
        
        title(sprintf('$N=%d$', N),'interpreter','latex');
        axis square;
        
        % Force axis boundary
        box on;
        
        filename = sprintf('../figs/derivative_N%d_order%d.pdf', N, order);
        print(f1, filename, '-dpdf');
        system(sprintf('pdfcrop --margins 10 %s %s', filename, filename));
        %close all;
        %}
        %% Find zeros
        %{
        assume(u,'positive');
        str2 = 'sol = solve(subs(hat(end), [';
        for i = 1:2^N-1
            str2 = [str2, sprintf('p%d ', i)];
        end
        str2 = [str2, '], [pval]) == 0, u)'];
        eval(str2)
        %}
    end
end

