function ChimeraDemoFast
% Produces metastable chimera-like states in a modular network of
% oscillators. The model is inspired by Abrams, et al., PRL 2008, and is
% parameterised by b and A, according to their paper.
% For more details of this system, see:
%
% Shanahan, M. (2010). Metastable Chimera States in Community-Structured
% Oscillator Networks. Chaos 20, 013108.
%
% Written by Murray Shanahan, August 2009 / December 2014

b = 0.1; % was 0.1
a = pi/2-b; % phase lag (alpha)

A = 0.2; % was 0.2
k1 = (1-A)/2; % inter-community coupling strength
k0 = 1-k1; % intra-community coupling strength

n0 = 32; n1 = 8; % community sizes

d0 = 32; d1 = 32; % numbers of connections at different community levels

N = n0*n1; % total number of oscillators
M = n1; % number of lowest level communities

% Build coupling matrix
K = zeros(N,N);
for i=1:N
   x1 = mod(ceil(i/n0)-1,n1)+1; % community number
   for j=i:N
      if i~=j % ignore diagonals
         y1 = mod(ceil(j/n0)-1,n1)+1; % community number
         if x1 == y1 % same community
            p = d0/n0;
            k = k0;
         else % different communities
            p = d1/(n0*n1);
            k = k1;
         end
         if rand < p
            K(i,j) = k;
            K(j,i) = k;
         end
      end
   end
end


w = ones(1,N); % identical natural frequencies


h = 0.05; % Runge-Kutta method step size

T = 800; % number of time steps

ws = 2; % window size for downsampling synchrony data


theta(1,:) = rand(1,N)*2*pi-pi; % random initial phases

Dmean = d0+d1; % average connections per oscillator


phi = zeros(T,M); % synchrony data

tot_phi = zeros(1,M); % running total for synchrony over window

Sync = zeros(T/ws,M); % downsampled synchrony data

hby2 = 0.5*h;

for t=2:T

   % Simulate the Kuramoto model
   temp1 = theta(t-1,:);
   for j=1:1/h
      temp2 = temp1;
      
      for i=1:N
         % Numerical simulation using 4th-order Runge-Kutta method
         rk1 = w(i);
         for k=1:N
            if K(i,k) ~= 0
               rk1 = rk1 + K(i,k)*sin(temp2(k)-temp1(i)-a)/Dmean;
            end
         end
         rk2 = w(i);
         for k=1:N
            if K(i,k) ~= 0
               rk2 = rk2 + K(i,k)*sin(temp2(k)-(temp1(i)+hby2*rk1)-a)/Dmean;
            end
         end
         rk3 = w(i);
         for k=1:N
            if K(i,k) ~= 0
               rk3 = rk3 + K(i,k).*sin(temp2(k)-(temp1(i)+hby2*rk2)-a)/Dmean;
            end
         end
         rk4 = w(i);
         for k=1:N
            if K(i,k) ~= 0
               rk4 = rk4 + K(i,k).*sin(temp2(k)-(temp1(i)+h*rk3)-a)/Dmean;
            end
         end
         temp1(i) = temp1(i) + h*(rk1+2*rk2+2*rk3+rk4)/6;
      end
            
   end
   theta(t,:) = temp1;

   % Compute synchrony within communities
   for i = 1:M
      for j = 1:n0
         x1 = theta(t,(i-1)*n0+j);
         phi(t,i) = phi(t,i)+exp(x1*sqrt(-1));
      end
   end
   phi(t,:) = abs(phi(t,:)/n0);
   
   theta(t,:) = mod(theta(t,:)+pi,2*pi)-pi; % normalise phases
   
   if mod(t,ws) == 0
      display(t)
      tot_phi = tot_phi + phi(t,:);
      Sync(t/ws,:) = tot_phi/ws; % average synchrony
      tot_phi = zeros(1,M);
      % Plot synchrony history
      subplot(3,M,2*M+1:3*M)
      plot(ws*(1:T/ws),Sync)
      xlabel('Time')
      ylabel('Synchrony')
      ylim([0 1])
   else
      tot_phi = tot_phi + phi(t,:);
   end

   % Plot phases
   Colours = lines(M);
   for i = 1:M
      subplot(3,M,[i i+M])
      plot(exp(theta(t,(i-1)*n0+1:i*n0)*sqrt(-1)),'.',...
         'Color',Colours(i,:),'MarkerSize',16)
      axis square
      xlim([-1 1])
      ylim([-1 1])
   end
   
   drawnow

end
