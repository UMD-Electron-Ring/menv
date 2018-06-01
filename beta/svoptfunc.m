function F =  svoptfunc( X )
% wrapper for optfunc; returns single valued min. function
%
    f = optfunc( X );
    F = f*f';
end