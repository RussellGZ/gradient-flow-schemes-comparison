clear; close all; rng(0);

% Parameters
eps2 = 0.001; h = 0.01; x = 0:h:1; N = numel(x);

% Finite-difference Laplacian with Neumann boundary conditions
Dh = zeros(N);
for i = 1:N
    if i == 1
        Dh(1,1) = -2; Dh(1,2) = 2;
    elseif i == N
        Dh(N,N-1) = 2; Dh(N,N) = -2;
    else
        Dh(i,i-1) = 1; Dh(i,i) = -2; Dh(i,i+1) = 1;
    end
end
Dh = Dh / h^2;

% Discrete energy
Eh = @(u) h * (-0.5*eps2*(u'*(Dh*u)) + sum((u.^2-1).^2)/4);

% Initial condition
u0 = 0.9*rand(N,1) + 0.05;
t_end = 40;

%% Scheme (2.6): maximum norm evolution
taus = [0.5,0.75,1,3];
figure;
for k = 1:4
    tau = taus(k);
    A   = eye(N) - tau*eps2*Dh;
    U   = u0;
    ts  = 0:tau:t_end;
    maxU = zeros(size(ts));
    for n = 1:length(ts)
        maxU(n) = max(abs(U));
        rhs = U + tau*(U - U.^3);
        U   = A \ rhs;
    end
    subplot(2,2,k);
    plot(ts, maxU, '-o'); grid on;
    xlabel('t'); ylabel('|\phi|_{max}');
    title(sprintf('$\\tau=%.2f$',tau),'Interpreter','latex');
end
sgtitle('Scheme (2.6)','Interpreter','latex');

%% Scheme (3.1): beta = 1
beta = 1; taus = [0.5,1,3];
figure;
for k = 1:4
    subplot(2,2,k);
    if k == 1
        hold on;
        for j = 1:length(taus)
            tau = taus(j);
            A   = (1 + beta*tau)*eye(N) - tau*eps2*Dh;
            U   = u0;
            ts  = 0:tau:t_end;
            Ehn = zeros(size(ts));
            for n = 1:length(ts)
                Ehn(n) = Eh(U);
                rhs    = (1+beta*tau)*U - tau*(U.^3 - U);
                U      = A \ rhs;
            end
            plot(ts, Ehn,'DisplayName',sprintf('\\tau=%.2f',tau));
        end
        grid on; ylabel('Energy');
        title(sprintf('Scheme (3.1), $\\beta=%d$: Energy',beta),'Interpreter','latex');
        legend('Location','northeast');
    else
        tau = taus(k-1);
        A   = (1 + beta*tau)*eye(N) - tau*eps2*Dh;
        U   = u0;
        ts  = 0:tau:t_end;
        maxU = zeros(size(ts));
        for n = 1:length(ts)
            maxU(n) = max(abs(U));
            rhs     = (1+beta*tau)*U - tau*(U.^3 - U);
            U       = A \ rhs;
        end
        plot(ts, maxU, '-o'); grid on;
        xlabel('t'); ylabel('|\phi|_{max}');
        title(sprintf('$\\tau=%.2f$',tau),'Interpreter','latex');
    end
end
sgtitle('Scheme (3.2), $\beta=1$','Interpreter','latex');

%% Scheme (3.1): beta = 2
beta = 2; % reuse taus
figure;
for k = 1:4
    subplot(2,2,k);
    if k == 1
        hold on;
        for j = 1:length(taus)
            tau = taus(j);
            A   = (1 + beta*tau)*eye(N) - tau*eps2*Dh;
            U   = u0;
            ts  = 0:tau:t_end;
            Ehn = zeros(size(ts));
            for n = 1:length(ts)
                Ehn(n) = Eh(U);
                rhs    = (1+beta*tau)*U - tau*(U.^3 - U);
                U      = A \ rhs;
            end
            plot(ts, Ehn,'DisplayName',sprintf('\\tau=%.2f',tau));
        end
        grid on; ylabel('Energy');
        title(sprintf('Scheme (3.1), $\\beta=%d$: Energy',beta),'Interpreter','latex');
        legend('Location','northeast');
    else
        tau = taus(k-1);
        A   = (1 + beta*tau)*eye(N) - tau*eps2*Dh;
        U   = u0;
        ts  = 0:tau:t_end;
        maxU = zeros(size(ts));
        for n = 1:length(ts)
            maxU(n) = max(abs(U));
            rhs     = (1+beta*tau)*U - tau*(U.^3 - U);
            U       = A \ rhs;
        end
        plot(ts, maxU, '-o'); grid on;
        xlabel('t'); ylabel('|\phi|_{max}');
        title(sprintf('$\\tau=%.2f$',tau),'Interpreter','latex');
    end
end
sgtitle('Scheme (3.2), $\beta=2$','Interpreter','latex');




function tight_pdf_export(figHandle, filename)
    set(figHandle, 'Units', 'Inches');
    pos = get(figHandle, 'Position');
    set(figHandle, 'PaperUnits', 'Inches');
    set(figHandle, 'PaperSize', [pos(3), pos(4)]);
    set(figHandle, 'PaperPosition', [0, 0, pos(3), pos(4)]);
    print(figHandle, filename, '-dpdf');
end

