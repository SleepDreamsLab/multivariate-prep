function count = poissonDraw(lambda)
% POISSONDRAW  One Poisson draw (Knuth's multiplication method).
% Adequate for the small rates used here (events per 30-s epoch) and avoids a dependency on
% the Statistics and Machine Learning Toolbox for poissrnd.
count = 0;
if lambda <= 0
    return
end
limit = exp(-lambda);
product = rand();
while product > limit
    count = count + 1;
    product = product * rand();
end
end
