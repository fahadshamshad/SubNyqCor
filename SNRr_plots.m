clear all; clc;
addpath('l1benchmark')
addpath('l1benchmark/l1benchmark')
addpath('phi2')
addpath('l1benchmark/l1benchmark/L1Solvers')
addpath('./yall1_rice/');
addpath('./sparsa/');


%parameters
SNRr = [5:5:40];
s = 10;
r = 10;%[18:4:52];
%measurements = [5:3:200];%[90:5:195];
RMSE_success = zeros(1,length(SNRr));

threshold_for_success = 10^-3;
M = 100; W = 1024; M1 = 20; %for rank
%R1 = 15;  %for column space
M2 = M-M1; %for sparsity
R2 = 60
R1 = floor(2.9*s*log10(W/s));

num_iter = 1;
error = zeros(1,length(num_iter));
rmse_error = zeros(1,length(SNRr));



for ss = 1:length(SNRr)
    signal_noise = SNRr(ss);
        for iter = 1:num_iter
             % for sparsity
            data = data_gene(s,r,M,W);

            %%%%
            tol=1e-6;
            A = AVMM(data,tol);
            A1 = A(1:M1,:);
            A2 = A(M1+1:end,:);

            %Column Space measurements
            D = binary(W);
            F = dftmtx(W)/sqrt(W);
            Q1 = accum(R1,W)*D*F;
            %Q1 = Q1;
            yc = data*Q1';
            [yc,sigma_c] = SNR_Set(yc,signal_noise);
            %noise_c = randn(size(yc));
            %noise_c = noise_c/norm(noise_c);
            %yc = yc+noise_c;

            %row space measurements
            Q2 = accum(R2,W)*D*F;
            Q2 = Q2;

            y = A2*data;
            f = Q2*y'; 
            [f,sigma_r] = SNR_Set(f,signal_noise);
            %noise = randn(size(f));
            %noise = noise/norm(noise);
            %f = f + noise;
            
            %y        = Phi*s; % noiseless measurements

            %----------------
            %%%%%%%% Measurment Matrix %%%%%%%%
            error = zeros(size(1,1));
            %%%%%%%%%%%%%%recovery L1 algorithm %%%%%%%%%%%;

            tic    
            s_r = sparse_recovery_yall1(Q2,f,W,M2,sigma_r);
            toc
            t=toc
            %%%% SVLS %%%% 

            [Xf, est_r] = SVLS(s_r',yc,A2,Q2',r);

            error_overall = RRMSE(Xf, data)
            error_sparse = RRMSE(s_r',y)
            error(iter) = error_overall; 
        end
       rmse_error(ss) =  mean(error);
end
rmse_error_dB = 10*log10(rmse_error.^2);
    
