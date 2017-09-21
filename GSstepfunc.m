function F =  GSstepfunc( X )
    f =  stepfunc2( X );
    F = rms(f);
end