function [SymExpOut,SubExpArray] = OptimiseAlgebra(SymExpIn,SubExpName)

% Loop through symbolic expression, identifying repeated expressions and
% bringing them out as shared expression or sub expressions
% do this until no further repeated expressions found
% This can significantly reduce computations

syms SubExpIn SubExpArray;

SubExpArray(1,1) = 'invalid';

index = 0;

f_complete = 0;

while f_complete == 0
    
    index = index + 1;

    % change variable naming scheme to avoid recent Matlab warnings
    SubExpIn = [SubExpName,'_l_',num2str(index),'_r_'];
    
    SubExpInStore{index} = SubExpIn;
    
    [SymExpOut,SubExpOut] = subexpr(SymExpIn,SubExpIn);

    % Fix crash when no sub-expression can be found
    if isempty(SubExpOut)
        
        break;
    end
    
    for k = 1:index
        
        if SubExpOut == SubExpInStore{k}
            f_complete = 1;
        end
        
    end
    
    if f_complete || index > 100
        
        SymExpOut = SymExpIn;
        
    else
        
        SubExpArray(index,1) = SubExpOut;
        
        SymExpIn = SymExpOut;
        
    end
end