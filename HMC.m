function [theta_hmc, acc, pot_energy] = HMC(epsilon_hmc, L_hmc, theta, A, Bpred, tBpred, Lambda, eta, sig, sum_B)

 q = theta;
 N = size(theta, 1);
 nbasis = size(theta, 2);

 % Independent standard normal variates 
 kinetic = normrnd(0, 1, N, nbasis);
 current_kinetic = kinetic;
 
 % Pre-compute matrices to speed up computation
 Leta = (Lambda * eta);
 AtBpred = A * tBpred;
		
 % Make a half step momentum at the beginning
 add_grad_q = - diag(sig) * Leta - sum_B';	
 ttheta = theta';
 grad_q = (AtBpred * exp(Bpred * ttheta)) +  diag(sig) * ttheta + add_grad_q;
 kinetic = kinetic - epsilon_hmc * (0.5)*(grad_q');
 
 % Time to compute (AtBpred * exp(Bpred * ttheta)) as a
 % function of N and p?
 % N =              100    250    500    750   1000    1250
 % p = 177       0.694s  1.302  2.279  3.557  4.633   6.830
 % p = 265       0.797s  1.727  3.001  4.548  6.429   7.914
 % p = 344       0.887s  1.855  3.227  4.724  6.448   8.257
% Bpred2 = dlmread('B.pred_less.txt');
% AtBpred2 = (Bpred2)';
% nbasis = size(Bpred2, 2);
% tic;
% testtheta = unifrnd(0, 1, [nbasis, 1000]);
% out = (AtBpred * exp(Bpred * testtheta));
% toc;
% xx = [100, 250, 500, 750, 1000, 1250];
% yy = [0.887  1.855  3.227  4.724  6.448   8.257];
% plot(xx, yy)
 
 % Alternate full steps for position and momentum
 for e = 1:L_hmc
			
    % Make a full step for the position
	q = q + epsilon_hmc * kinetic;
    tq = q';
			
	% Make a full step for the momentum, except at end of trajectory
	if e < L_hmc
		grad_q = (AtBpred * exp(Bpred * tq)) + diag(sig) * (tq) + add_grad_q;
		kinetic = kinetic - epsilon_hmc * grad_q';
    end
 end
		
 % Make a half step momentum at the end
 grad_q = (AtBpred * exp(Bpred * (tq))) + diag(sig) * (tq) + add_grad_q;
 kinetic = kinetic - epsilon_hmc * 0.5*(grad_q');
	
 % Negate momentum at end of trajectory to make the proposal symmetric
 kinetic = - kinetic;
		
 % Evaluate potential and kinetic energies at start and end of trajectory
 current_U = A.*sum(exp(Bpred * ttheta))' - diag(sum_B * ttheta) + 0.5 .* diag((theta - (Leta)') * diag(sig) * (theta - (Leta)')');
 current_K = sum((current_kinetic.^2)')'/2;
 proposed_U = A.*sum(exp(Bpred * tq))' - diag(sum_B * tq) + 0.5 .* diag((q - (Leta)') * diag(sig) * (q - (Leta)')');
 proposed_K = sum((kinetic.^2)')'/2;
		
 % Accept or reject the state at end of trajectory, returning either the position 
 % at the end of the trajectory or the initial position
		
 pratio = current_U - proposed_U + current_K - proposed_K;	% Log acceptance ratio
 u = log(rand([N, 1]));
 acc = 1 * (u < pratio);
 ind = find(acc == 1);
 theta_hmc = theta; theta_hmc(ind, :) = q(ind, :); % theta_hmc = old, excpet for those proposal which get accepted 
		
 pot_energy = A*sum(exp(Bpred * (theta_hmc)'), 1)' - diag(sum_B * (theta_hmc)') + ... 
    0.5 * diag((theta_hmc - (Leta)') * diag(sig) * (theta_hmc - (Leta)')');
end