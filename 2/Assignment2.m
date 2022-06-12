clear all;
close all;

%% Variables
E = 1;
s = randsrc(10^5, 1);
SNR_db =(-6 : 1 : 6);
SNR = 10.^(SNR_db / 20);
N = E./SNR;
var = sqrt(N / 2);
norm = 10^5;


%% Gaussian noise
SG = zeros(10^5, 13);
for i = 1 : 13
   WG(:, i) = var(i) * randn(10^5, 1); 
   RG(:, i) = s + WG(:, i);
   for j = 1 : 10^5   
       if (RG(j, i) > 0)
           SG(j, i) = 1;
       else
           SG(j, i) = -1;
       end
   end
end

prob_error_A = zeros(13,1);

for i = 1 : 13
   for j = 1 : 10^5   
       if(SG(j, i) - s(j) ~= 0)
          prob_error_A(i) = prob_error_A(i) + 1; 
       end  
   end
end

prob_error_gauss = prob_error_A / norm;

%Graph
figure(1);
plot(SNR_db, prob_error_gauss / norm ,'-h');
grid on
title('P_e_r_r_o_r - Gaussian noise');
xlabel('SNR[dB]');
ylabel('P_e(SNR[dB])');
movegui('west');


%% Laplacian noise
SL=zeros(10^5,13);
for i = 1 : 13
   WL(:, i) = laprnd(10^5, 1, 0, var(i)); 
   RL(:, i) = s + WL(:, i);
   for j = 1 : 10^5   
       if (RL(j, i) > 0)
           SL(j, i) = 1;
       else
           SL(j, i) = -1;
       end
   end
end

prob_error_laplace = zeros(13, 1);

for i = 1 : 13
    for j = 1 : 10^5   
        if(SL(j, i) - s(j) ~= 0)
           prob_error_laplace(i) = prob_error_laplace(i) + 1; 
        end  
    end
end

%Graph:
figure(2);
plot(SNR_db, prob_error_laplace/norm, 'r-h');
grid on
title('P_e_r_r_o_r - Laplacian noise');
xlabel('SNR[dB]');
ylabel('P_e(SNR[dB])');
movegui('east');

%% Compare Plots
figure(3);
plot(SNR_db, prob_error_A/norm, 'b-*');
hold on
plot(SNR_db, prob_error_laplace/norm, 'r-o');
grid on
title('P_e_r_r_o_r - Comparison');
xlabel('SNR[dB]');
ylabel('P_e(SNRdB)');
legend('Gaussian noise','Laplacian noise')
movegui('north');

%% Sub-Functions
function y  = laprnd(m, n, mu, sigma)
%LAPRND generate i.i.d. laplacian random number drawn from laplacian distribution
%   with mean mu and standard deviation sigma. 
%   mu      : mean
%   sigma   : standard deviation
%   [m, n]  : the dimension of y.
%   Default mu = 0, sigma = 1. 
%   For more information, refer to
%   http://en.wikipedia.org./wiki/Laplace_distribution

%   Author  : Elvis Chen (bee33@sjtu.edu.cn)
%   Date    : 01/19/07

%Check inputs
if nargin < 2
    error('At least two inputs are required');
end

if nargin == 2
    mu = 0; sigma = 1;
end

if nargin == 3
    sigma = 1;
end

% Generate Laplacian noise
u = rand(m, n)-0.5;
b = sigma / sqrt(2);
y = mu - b * sign(u).* log(1- 2* abs(u));
end





