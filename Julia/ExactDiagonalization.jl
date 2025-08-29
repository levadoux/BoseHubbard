using LinearAlgebra
using SparseArrays
using Combinatorics

function size_Fock_space(N, M)
    return binomial(N + M - 1, N)
end

function state_to_positions(N, M, state)
    positions = Int[]
    for j in 1:M
        append!(positions, fill(M - j + 1, state[M - j + 1]))
    end
    return positions
end

function positions_to_label(N, M, positions)
    label = 1
    for j in 1:N
        label += size_Fock_space(j, M - positions[j])
    end
    return label
end

function labeling(N, M, state)
    positions = state_to_positions(N, M, state)
    return positions_to_label(N, M, positions)
end

function label_to_positions(N, M, label)
    positions = zeros(Int, N)
    rest = label
    for j in 0:N-1
        exp = 0
        while size_Fock_space(N - j, exp) < rest
            exp += 1
        end
        m = M + 1 - exp  # car exp = M - m + 1 à la fin
        positions[N - j] = m
        rest -= size_Fock_space(N - j, M - m)
    end
    return positions
end

function positions_to_state(N, M, positions)
    state = zeros(Int, M)
    for site in positions
        state[site] += 1
    end
    return state
end

function unlabeling(N, M, label)
    positions = label_to_positions(N, M, label)
    return positions_to_state(N, M, positions)
end

function creation_annihilation(i, j, state)
    new_state = copy(state)
    n_i = new_state[i]
    n_j = new_state[j]
    if i == j
        return new_state, new_state[i]
    elseif n_j > 0
        new_state[i] += 1
        new_state[j] -= 1
        coeff = sqrt(n_j) * sqrt(n_i + 1)
        return new_state, coeff
    else
        return state, 0.0
    end
end

function BoseHubbard_Hj(N, M)
    dim = size_Fock_space(N, M)
    Hj = zeros(Float64, dim, dim)
    for label in 1:dim
        state = unlabeling(N, M, label)
        for j in 1:M-1
            new_state1, coeff1 = creation_annihilation(j, j+1, state)
            new_state2, coeff2 = creation_annihilation(j+1, j, state)
            if coeff1 != 0
                new_label1 = labeling(N, M, new_state1)
                Hj[new_label1, label] = coeff1
            end
            if coeff2 != 0
                new_label2 = labeling(N, M, new_state2)
                Hj[new_label2, label] = coeff2
            end
        end
        if M > 2
            new_state1, coeff1 = creation_annihilation(M, 1, state)
            new_state2, coeff2 = creation_annihilation(1, M, state)
            if coeff1 != 0
                new_label1 = labeling(N, M, new_state1)
                Hj[new_label1, label] = coeff1
            end
            if coeff2 != 0
                new_label2 = labeling(N, M, new_state2)
                Hj[new_label2, label] = coeff2
            end
        end
    end
    return Hj
end

function BoseHubbard_Hu(N, M)
    dim = size_Fock_space(N, M)
    Hu = zeros(Float64, dim, dim)
    for label in 1:dim
        state = unlabeling(N, M, label)
        Hu[label, label] = sum(0.5 * n * (n - 1) for n in state)
    end
    return Hu
end

function Single_Hopping(N, M, i, j)
    dim = size_Fock_space(N, M)
    H_ij = zeros(Float64, dim, dim)
    for label in 1:dim
        state = unlabeling(N, M, label)
        new_state, coeff = creation_annihilation(i, j, state)
        if coeff != 0
            new_label = labeling(N, M, new_state)
            H_ij[new_label, label] = coeff
        end
    end
    return H_ij
end

function OBDM_ground(N, M, U, J)
    dim = size_Fock_space(N, M)
    H = -J * BoseHubbard_Hj(N, M) + U * BoseHubbard_Hu(N, M)
    evals, evecs = eigen(H)
    ground = evecs[:, argmin(evals)]
    rho1 = zeros(Float64, M, M)
    for i in 1:M, j in 1:M
        H_ij = Single_Hopping(N, M, i, j)
        rho1[i, j] = ground' * H_ij * ground
    end
    return rho1
end

function crea_crea_ann_ann(k, l, j, i, state)
    new_state = copy(state)
    n_i = new_state[i]
    n_j = new_state[j]
    if i == j && n_i > 1
        new_state[i] -= 2
        n_l = new_state[l]
        new_state[l] += 1
        n_k = new_state[k]
        new_state[k] += 1
        coeff = sqrt(n_i * (n_i - 1)) * sqrt((n_l + 1) * (n_k + 1))
    elseif i != j && n_i > 0 && n_j > 0
        new_state[i] -= 1
        new_state[j] -= 1
        n_l = new_state[l]
        new_state[l] += 1
        n_k = new_state[k]
        new_state[k] += 1
        coeff = sqrt(n_i * n_j) * sqrt((n_l + 1) * (n_k + 1))
    else
        return state, 0.0
    end
    return new_state, coeff
end

function Pair_Hopping(k, l, j, i, N, M)
    dim = size_Fock_space(N, M)
    Hijkl = zeros(Float64, dim, dim)
    for label in 1:dim
        state = unlabeling(N, M, label)
        new_state, coeff = crea_crea_ann_ann(k, l, j, i, state)
        if coeff != 0
            new_label = labeling(N, M, new_state)
            Hijkl[new_label, label] = coeff
        end
    end
    return Hijkl
end

function TBDM_ground(N, M, U, J)
    dim = size_Fock_space(N, M)
    H = -J * BoseHubbard_Hj(N, M) + U * BoseHubbard_Hu(N, M)
    evals, evecs = eigen(H)
    ground = evecs[:, argmin(evals)]
    rho2 = zeros(Float64, M*M, M*M)
    for i in 1:M, j in 1:M, k in 1:M, l in 1:M
        Hijkl = Pair_Hopping(k, l, j, i, N, M)
        rho2[(i-1)*M+j, (k-1)*M+l] = ground' * Hijkl * ground
    end
    return rho2
end

