using JuMP, SCS, LinearAlgebra

function groundSDP(N, M, U, J)
    # Construct the one-body Hamiltonian H1 (M x M)
    H1 = zeros(M, M)
    for i in 1:(M-1)
        H1[i, i+1] = -J   # upper diagonal
        H1[i+1, i] = -J   # lower diagonal
    end
    # Periodic boundary conditions
    H1[1, M] = -J
    H1[M, 1] = -J

    # Construct the two-body Hamiltonian H2 (M^2 x M^2)
    H2 = zeros(M*M, M*M)
    for i in 1:M
        idx = (i-1)*M+i 
        H2[idx, idx] = U/2
    end

    # Create the optimization model with the SCS solver
    model = Model(SCS.Optimizer)
    set_optimizer_attribute(model, "verbose", 0)

    # Variables:
    # Gamma1 is a real symmetric PSD matrix of size (M x M)
    @variable(model, Gamma1[1:M, 1:M], PSD)
    # Gamma2 is a real symmetric PSD matrix of size (M^2 x M^2)
    @variable(model, Gamma2[1:M*M, 1:M*M], PSD)
    # Define G as a PSD matrix (M^2 x M^2)
    @variable(model, G[1:M*M, 1:M*M], PSD)

    # Objective: minimize Tr(H1 * Gamma1) + Tr(H2 * Gamma2)
    obj = sum(H1[i,j] * Gamma1[j,i] for i in 1:M, j in 1:M) +
          sum(H2[i,j] * Gamma2[j,i] for i in 1:(M*M), j in 1:(M*M))
    @objective(model, Min, obj)

    # Constraints:
    # Trace(Gamma1) = N
    @constraint(model, sum(Gamma1[i,i] for i in 1:M) == N)
    # Trace(Gamma2) = N(N - 1)
    @constraint(model, sum(Gamma2[i,i] for i in 1:(M*M)) == N * (N - 1))

    # Contraction constraint:
    # For all i, j in {1,...,M}, sum_k Gamma2[(i-1)*M + k, (j-1)*M + k] == (N-1) * Gamma1[i,j]
    for i in 1:M
        for j in 1:M
            @constraint(model, sum(Gamma2[(i-1)*M + k, (j-1)*M + k] for k in 1:M) == (N - 1) * Gamma1[i,j])
        end
    end

    # Symmetry constraint on Gamma2:
    # For all i, j, k in {1,...,M}, and l in {1,...,k-1}:
    # Gamma2[(i-1)*M + j, (k-1)*M + l] == Gamma2[(i-1)*M + j, (l-1)*M + k]
    for i in 1:M
        for j in 1:M
            for k in 1:M
                for l in 1:(k-1)
                    @constraint(model, Gamma2[(i-1)*M + j, (k-1)*M + l] == Gamma2[(i-1)*M + j, (l-1)*M + k])
                end
            end
        end
    end    

    # Constraint linking G, Gamma1, and Gamma2:
    # G[(i-1)*M + j, (k-1)*M + l] = δ_{jl} * Gamma1[i,k] + Gamma2[(i-1)*M + l, (k-1)*M + j]
    for i in 1:M
        for j in 1:M
            for k in 1:M
                for l in 1:M
                    delta_jl = (j == l) ? 1.0 : 0.0
                        @constraint(model, G[(i-1)*M + j, (k-1)*M + l] == delta_jl * Gamma1[i, k] + Gamma2[(i-1)*M + l, (k-1)*M + j])
                end
            end
        end
    end

    # Solve the SDP
    optimize!(model)

    # Return optimal values of Gamma1 and Gamma2
    return value.(Gamma1), value.(Gamma2), value.(obj)
end
