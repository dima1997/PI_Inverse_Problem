function [noise_power]=discrepancy_principle(Noise,lambda_range)
    noise_power=norm(Noise,'fro')^2*ones(length(lambda_range),1);
end