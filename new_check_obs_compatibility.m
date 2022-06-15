function [compatible,tF] = new_check_obs_compatibility(ftruth,fobs)
% Check if truth and obs are compatible
compatible = false;
tF = 0.0;

% Truth file e.g. L04M2_N960_K32_F8.00_I1_b1.00_c1.00_tf0.05_spinup100_tsteps10000_seed51422
% Obs file e.g. obs_tf0.05_nc10000_R1.00_N960_K32_F8.00_I1_b1.00_c1.00_tseed51422_oseed73033
struth=strsplit(ftruth,'_');
sobs=strsplit(fobs,'_');

% Compare time steps
tf_truth = struth{9};
tf_obs = sobs{2};
if ~strcmp(tf_truth,tf_obs)
    return
end
tF = str2num(tf_obs(3:end));
% Compare seeds
seed_truth = ['t' struth{12}(1:end-4)];
seed_obs = sobs{11};
if ~strcmp(seed_truth,seed_obs)
    return
end
% Compare grid size
ng_truth = struth{3}(1:end);
ng_obs = sobs{5}(2:end);
if ng_obs ~= ng_truth
    return
end
% Compare number of cycles
ncycles_truth = struth{11}(7:end);
ncycles_obs = sobs{3}(3:end);
if ncycles_obs > ncycles_truth
    return
elseif ncycles_obs < ncycles_truth
    fprintf('Warning: ncycles_obs = %d < ncycles_truth = %d\n',...
        ncycles_obs,ncycles_truth);
end

% All checks passed
compatible = true;
