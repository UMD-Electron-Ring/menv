function F =  GSstepfunc( X )
    f =  stepfunc2( X );
    F = sum(abs(f));
end