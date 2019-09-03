function y = logTrafo(x, direction)

if(~exist('direction','var') || isempty(direction))
    direction = 1;
end

if direction==1
    if(istable(x))
        y = array2table(log2(1+x{:,:}));
        y.Properties.RowNames = x.Properties.RowNames;
        y.Properties.VariableNames = x.Properties.VariableNames;
        clear x
    else
        y = log2(1+x);
    end
elseif direction==-1
    if(istable(x))
        y = array2table((2.^x{:,:})-1);
        y.Properties.RowNames = x.Properties.RowNames;
        y.Properties.VariableNames = x.Properties.VariableNames;
        clear x
    else
        y = (2.^x)-1;
    end
end
