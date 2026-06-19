function digits = killN1(scoringDigits)
% Replace isolated N1 epochs at stage boundaries with their neighbour's stage.
    digits = scoringDigits;
    while any(digits == -1)
        firstN1 = find([0,  diff(digits == -1)] ==  1);
        lastN1  = find([diff(digits == -1),  0] == -1);
        digits(firstN1) = digits(firstN1 - 1);
        digits(lastN1)  = digits(lastN1  + 1);
    end
end
