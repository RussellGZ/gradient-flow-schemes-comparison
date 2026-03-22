clear; close all; rng(0);
L = 2*pi;
N = 512; 
dx = L/N;
x = (0:N-1)*dx;
[X,Y] = meshgrid(x,x);
epsilon = 0.1; 
dt = 1e-3; 
Tend = 8;
M = round(Tend/dt);

k = [0:N/2-1 -N/2:-1] * (2*pi/L);
[KX,KY]= meshgrid(k,k);
K2 = KX.^2 + KY.^2;
K4 = K2.^2;


Lspec = -epsilon^2 * K4;
expL = exp(dt * Lspec);
phi1 = (expL - 1) ./ Lspec;
phi2 = (expL - 1 - dt*Lspec) ./ (Lspec.^2);
phi1(1,1) = dt;
phi2(1,1) = dt^2/2;

u = 0.05 * sin(X) .* sin(Y);
snapshots = zeros(N,N,4);
snap_times = [2,4,6,8];
snap_idx = round(snap_times/dt);

E = zeros(1, M+1);
E(1) = compute_energy(u, KX, KY, dx, epsilon);
for n = 1:M
    NL0_hat = -K2 .* fft2(u.^3 - u);
    u_hat = fft2(u);
    v_hat = expL .* u_hat + phi1 .* NL0_hat;
    v = real(ifft2(v_hat));
    NL1_hat = -K2 .* fft2(v.^3 - v);
    u_hat_new = expL .* u_hat + phi1 .* NL0_hat + phi2 .* (NL1_hat - NL0_hat);
    u = real(ifft2(u_hat_new));
    if ismember(n, snap_idx)
        j = find(snap_idx==n,1);
        snapshots(:,:,j) = u;
    end
    E(n+1) = compute_energy(u, KX, KY, dx, epsilon);
end
figure('Position',[200 200 800 800]);
for j = 1:4
    subplot(2,2,j);
    imagesc(x, x, snapshots(:,:,j));
    axis equal tight off;
    title(sprintf('T = %d', snap_times(j)));
    colormap(parula);
end
figure;
plot(0:dt:Tend, E, 'LineWidth', 1.5);
xlabel('Time, t');
ylabel('Free energy E(t)');
title('Energy evolution for Cahn–Hilliard via ETDRK2');
grid on;
function E = compute_energy(u, KX, KY, dx, epsilon)
    Uhat = fft2(u);
    ux   = real(ifft2(1i * KX .* Uhat));
    uy   = real(ifft2(1i * KY .* Uhat));
    Ekin = 0.5 * sum((ux.^2 + uy.^2), 'all') * dx^2;
    Epot = sum((1 - u.^2).^2, 'all')/(4*epsilon^2) * dx^2;
    E = Ekin + Epot;
end
