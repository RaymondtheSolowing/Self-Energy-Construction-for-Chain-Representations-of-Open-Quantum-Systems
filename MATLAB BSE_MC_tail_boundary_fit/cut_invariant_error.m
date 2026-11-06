function chi=cut_invariant_error(SigmaA,SigmaB)
%CUT_INVARIANT_ERROR chi=|Sigma_a-Sigma_b|/|Sigma_a+conj(Sigma_b)|.
% Both self-energies are presented to the same site at the same p.
% On the real axis, chi is unchanged by every step of the self-energy
% recurrence, so its value at the cut equals its value at the system port.
% It is the magnitude of the reflection amplitude of the closure.
assert(isequal(size(SigmaA),size(SigmaB)),'Self-energies must have equal sizes.');
chi=abs(SigmaA-SigmaB)./abs(SigmaA+conj(SigmaB));
end
