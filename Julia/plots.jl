using LinearAlgebra
using Arpack
using Plots

# Display functions for exact diagonalization

function curve_ground_state_energy(N, M, alpha)
    g = Float64[]
    Hj = BoseHubbard_Hj(N, M)
    Hu = BoseHubbard_Hu(N, M)
    for U in alpha
        H = -Hj + U * Hu
        Eg = eigs(H; nev=1, which=:SR)[1]
        push!(g, real(Eg))
    end
    return g
end

function plot_ground_state_energy(N, M)
    alpha = range(0, stop=100, length=100)
    g = curve_ground_state_energy(N, M, alpha)
    plot(alpha, g, xlabel="\$U/J\$", ylabel="\$E_{gs}/J\$",
         title="Energy of the ground state for N = $N, M = $M", legend=false)
end

function plot_OBDM(N, M, U, J)
    rho1 = OBDM_ground(N, M, U, J)
    heatmap(rho1, color=:plasma, clim=(0, N/M),
            title="N = $N, M = $M, U = $U J",
            aspect_ratio=:equal)
end

function curve_condensed_fraction(N, M, alpha)
    condensed_fraction = Float64[]
    for U in alpha
        rho1 = OBDM_ground(N, M, U, 1)
        c = maximum(real(eigen(rho1).values))
        c = c / tr(rho1)
        push!(condensed_fraction, c)
    end
    return condensed_fraction
end

function plot_condensed_fraction(N, M)
    alpha = range(0, stop=100, length=100)
    cf = curve_condensed_fraction(N, M, alpha)
    plot(alpha, cf, xlabel="\$U/J\$", ylabel="\$f_{cond}\$",
         title="Condensed fraction for N = $N, M = $M", legend=false)
end

# Display functions for SDP

function curve_ground_state_energy_SDP(N, M, alpha::Vector{Float64})
    g_SDP = Float64[]
    for U in alpha
        try
            start = time()
            Gamma1, Gamma2, energy = groundSDP(N, M, U, 1.0)          
            push!(g_SDP, energy)
        catch err
            push!(g_SDP, 0.0)
            @warn "Failed for U = $U" exception=err
        end
    end
    return g_SDP
end

function plot_ground_state_energy_SDP(N, M)
    alpha = range(0, stop=100, length=100)
    g_SDP = curve_ground_state_energy_SDP(N, M, alpha)
    plot(alpha, g_SDP, xlabel="\$U/J\$", ylabel="\$E_{gs}/J\$",
         title="Energy of the ground state for N = $N, M = $M (SDP)", legend=false)
end

function plot_OBDM_SDP(N, M, U, J)
    Gamma1, Gamma2, energy = groundSDP(N, M, U, 1.0)
    rho1_SDP = Gamma1
    heatmap(rho1_SDP, color=:plasma, clim=(0, N/M),
            title="N = $N, M = $M, U = $U J (SDP)",
            aspect_ratio=:equal)
end

function curve_condensed_fraction_SDP(N, M, alpha)
    condensed_fraction = Float64[]
    for U in alpha
        Gamma1, Gamma2, energy = groundSDP(N, M, U, 1.0)
        rho1_SDP = Gamma1
        c = maximum(real(eigen(rho1_SDP).values))
        c = c / tr(rho1_SDP)
        push!(condensed_fraction, c)
    end
    return condensed_fraction
end

function plot_condensed_fraction_SDP(N, M)
    alpha = range(0, stop=100, length=100)
    cf = curve_condensed_fraction_SDP(N, M, alpha)
    plot(alpha, cf, xlabel="\$U/J\$", ylabel="\$f_{cond}\$",
         title="Condensed fraction for N = $N, M = $M (SDP)", legend=false)
end

# Display functions to compare SDP and exact diagonalization

function compare_ground_state_energy(N, M)
    alpha = range(0, stop=100, length=100)
    g_exact = curve_ground_state_energy(N, M, alpha)
    g_SDP   = curve_ground_state_energy_SDP(N, M, alpha)

    plot(alpha, g_exact, label="Exact", xlabel="\$U/J\$", ylabel="\$E_{gs}/J\$",
         title="Ground state energy for N = $N, M = $M")
    plot!(alpha, g_SDP, label="SDP")
end

function compare_condensed_fraction(N, M)
    alpha = range(0, stop=100, length=100)
    cf_exact = curve_condensed_fraction(N, M, alpha)
    cf_SDP   = curve_condensed_fraction_SDP(N, M, alpha)

    plot(alpha, cf_exact, label="Exact", xlabel="\$U/J\$", ylabel="\$f_{cond}\$",
         title="Condensed fraction for N = $N, M = $M")
    plot!(alpha, cf_SDP, label="SDP")
end

# To evaluate the accuracy of the SDP we use the trace distance between the normalized OBDMs.
# Computes the trace distance between two matrices A and B
# trace_distance(A,B) = 0.5 * ||A - B||_1, where ||.||_1 is the trace norm (sum of singular values)
function trace_distance(A, B)
    diff = A - B                      
    return 0.5 * sum(svdvals(diff))     # Sum of singular values gives the trace norm
end

# We define the accuracy by 1-trace_distance
function SDP_accuracy(N, M, U, J)
    rho1_normalized = OBDM_ground(N, M, U, J)
    rho1_SDP, _, _ = groundSDP(N, M, U, J)
    rho1_SDP_normalized = rho1_SDP * M / N
    return 1-trace_distance(rho1_normalized, rho1_SDP_normalized)
end
    


