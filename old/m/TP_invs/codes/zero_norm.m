function z = zero_norm(x)
    z = sum(x > 0.01*max(x));
end