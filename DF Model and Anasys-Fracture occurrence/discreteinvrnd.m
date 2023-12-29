function X = discreteinvrnd(p,m,n)

  X = zeros(m,n); % Preallocate memory
  for i = 1:m*n
    u = rand;
    I = find(u < cumsum(p));
    X(i) = min(I);
  end

end