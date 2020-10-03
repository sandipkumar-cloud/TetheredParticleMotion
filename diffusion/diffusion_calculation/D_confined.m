function [square_deviation] = D_confined(t,L,var_sigma)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
A=0;
n=0;
b=-1;
    while (A-b)> (2*eps)
        b=A;
        A=A+(1/((2*n+1)^4))*exp(-((t/2)*((((2*n+1)*pi*var_sigma)/L)^2)));
        n=n+1;
    end
    
   square_deviation=(L^2)/6-((16*(L^2))/(pi^4))*A;
end



