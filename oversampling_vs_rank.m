clear all; clc;
addpath('l1benchmark')
addpath('l1benchmark/l1benchmark')
addpath('phi2')
addpath('l1benchmark/l1benchmark/L1Solvers')
addpath('./yall1_rice/');


%parameters

s = 15;
rank = 5;%[18:4:52];
%R1 = 5;%[90:5:195];
R1_success = zeros(1,length(rank));

threshold_for_success = 10^-3;
M = 100; W = 1024; M2 = 20; %for rank
%R1 = 15;  %for column space
M1 = M-M2; %for sparsity
R2 = 80;

eta = 15;

R1 = in_over_out_r1(eta,s,rank,M,M1,M2,R2)

num_iter = 3;
for ss = 1:length(rank)
    r = rank(ss);
    
    for meas = 1:1
        num_success = 0;
        for iter = 1:num_iter
           
            %R1 = measurements(meas)%floor(2.9*s*log10(W/s)); % for sparsity

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

            %row space measurements
            Q2 = accum(R2,W)*D*F;
            Q2 = Q2;

            y = A2*data;
            f = Q2*y';

            %----------------
            %%%%%%%% Measurment Matrix %%%%%%%%
            error = zeros(size(1,1));
            %%%%%%%%%%%%%%recovery L1 algorithm %%%%%%%%%%%;

            tic    
            s_r = sparse_recovery_yall1(Q2,f,W,M2,0);
            toc
            t=toc
            %%%% SVLS %%%% 

            [Xf, est_r] = SVLS(s_r',yc,A2,Q2',r);

            error_overall = RRMSE(Xf, data)
            error_sparse = RRMSE(s_r',y)

            num_success = num_success + (error_overall < threshold_for_success);
            R1_success(ss) = R1;
            R1_success
            
        end
    end
end
    
