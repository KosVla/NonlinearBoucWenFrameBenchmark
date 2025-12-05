function kbw = compute_kbw(zeta, delta_u, Alpha, Beta, Gamma, N, ve, ne)
    kbw = (Alpha - ve.*(Beta.*sign(zeta).*sign(delta_u) + Gamma).* abs(zeta).^(N))./ne;
end