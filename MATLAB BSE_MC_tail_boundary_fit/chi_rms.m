function value=chi_rms(w,chi)
%CHI_RMS Root mean square of the cut-invariant error over a frequency grid.
value=sqrt(trapz(w(:),chi(:).^2)/(max(w)-min(w)));
end
