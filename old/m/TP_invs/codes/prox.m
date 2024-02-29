function y = prox(b, x)
    y = (x-b) .* (x > b) + (b+x) .* (x < -b);
end