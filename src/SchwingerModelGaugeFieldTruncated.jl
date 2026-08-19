##############################################################################################
# Codes to explore various ways of truncating the gauge degrees of freedom in the Schwinger
# model. The codes are based on ITensor and its DMRG code
##############################################################################################
using LinearAlgebra
using SparseArrays
using ITensors
using ITensorMPS
using Printf

"""
    ITensors.op(::OpName"L", ::SiteType"Qudit", dim::Int)

The electric field operator as an ITensor operator
"""
function ITensors.op(::OpName"L", ::SiteType"Qudit", dim::Int)
    entries = collect((-(dim-1)/2):1:((dim-1)/2))
    if iseven(dim)
        entries = collect((-dim/2):1:(dim/2))
    end
    return diagm(entries)
end

"""
    ITensors.op(::OpName"Lsq", ::SiteType"Qudit", dim::Int)

The electric field operator squared as an ITensor operator
"""
function ITensors.op(::OpName"Lsq", ::SiteType"Qudit", dim::Int)
    L = Matrix(op("L", siteind("Qudit", dim=dim)).tensor)
    return L * L
end

"""
    ITensors.op(::OpName"U", ::SiteType"Qudit", dim::Int)

The U operator as an ITensor operator
"""
function ITensors.op(::OpName"U", ::SiteType"Qudit", dim::Int)
    entries = ones(dim - 1)
    return diagm(-1 => entries)
end

"""
    ITensors.op(::OpName"Ud", ::SiteType"Qudit", dim::Int)

The hermitian conjugate of U as an ITensor operator
"""
function ITensors.op(::OpName"Ud", ::SiteType"Qudit", dim::Int)
    U = Matrix(op("U", siteind("Qudit", dim=dim)).tensor)
    return U'
end

"""
    ITensors.op(::OpName"Qodd", ::SiteType"S=1/2")

The staggered charge for odd matter sites
"""
function ITensors.op(::OpName"Qodd", ::SiteType"S=1/2")
    Z = Matrix(op("Z", siteind("S=1/2")).tensor)
    Q = 0.5 * (Z - I)
    return Q
end

"""
    ITensors.op(::OpName"Qeven", ::SiteType"S=1/2")

The staggered charge for even matter sites
"""
function ITensors.op(::OpName"Qeven", ::SiteType"S=1/2")
    Z = Matrix(op("Z", siteind("S=1/2")).tensor)
    Q = 0.5 * (Z + I)
    return Q
end

"""
    ITensors.op(::OpName"Qoddsq", ::SiteType"S=1/2")

The square of the staggered charge for odd matter sites
"""
function ITensors.op(::OpName"Qoddsq", ::SiteType"S=1/2")
    Q = Matrix(op("Qodd", siteind("S=1/2")).tensor)
    return Q * Q
end

"""
    ITensors.op(::OpName"Qevensq", ::SiteType"S=1/2")

The square of the staggered charge for even matter sites
"""
function ITensors.op(::OpName"Qevensq", ::SiteType"S=1/2")
    Q = Matrix(op("Qeven", siteind("S=1/2")).tensor)
    return Q * Q
end


"""
    getSites(N::Int, dim::Int)

Provide the definition of the lattice where dim specifies to which dimension the gauge degrees of freedom are truncated
"""
function getSites(N::Int, d::Int)
    Ntotal = 2 * N - 1
    s = Vector{Index{Int64}}(undef, Ntotal)
    # Make an array of 'site' indices alternating between spin 1/2 and gauge links with dimension dim
    for i = 1:Ntotal
        if isodd(i)
            s[i] = siteind("S=1/2", i)
        else
            s[i] = siteind("Qudit", i, dim=d)
        end
    end
    return s
end


############################################################################################
# The Hamiltonian
############################################################################################

"""
    getSchwingerHamiltonianMPO(J::Real, mu::Real, epsilon::Real, l0::Real, lambda::Real, eta::Real, s::Vector{Index{Int64}}, q::Vector{Int}=Vector{Int64}(undef, 0))

Given a lattice, build the Hamiltonian MPO for the Schwinger model. The coefficients are

``H = J H_{\\mathrm{kin}} + m H_{m} + epsilon H_{\\mathrm{el}} + \\lambda P_{\\mathrm{Gauss}} + \\eta P_{Q}``

where the individual terms are the kinetic term, the mass term, the electric energy term and a penalty enforcing gauge invariance and vanishing total charge.
The individual Hamiltonian terms correspond to

``H_{\\mathrm{kin}} = \\frac{1}{2a}\\sum_{l}\\left(\\phi_{l}^\\dagger U \\phi_{l+1} + ext{h.c.}\\right)``

``H_{m} = \\sum_{l}(-1)^l \\phi^\\dagger_{l}\\phi_{l}``

``H_{\\mathrm{el}}} = \\sum_l L_{l,l+1}^2``

To obtain the dimensionless Schwinger Hamiltonian the coefficients have to chosen as

``J = x, mu = 2 m/g \\times \\sqrt{x}, \\epsilon = 1``
"""
function getSchwingerHamiltonianMPO(J::Real, mu::Real, epsilon::Real, l0::Real, lambda::Real, eta::Real, s::Vector{Index{Int64}})
    # Compute the system size
    N = Int(round((length(s) + 1) / 2))
    @assert(2 * N - 1 == length(s))
    # Build the Hamiltonian MPO, seemingly unnecessary if-statements, are there to prevent terms in the MPO with prefactors that are numerical zeros and could potentially increase the bond dimension without having any effect            
    terms = OpSum()
    for n = 1:2:length(s)
        nsite = Int(round((n + 1) / 2))

        if n < length(s) - 1
            # The kinetic term
            if J != 0
                terms -= J, "S+", n, "U", n + 1, "S-", n + 2
                terms -= J, "S-", n, "Ud", n + 1, "S+", n + 2
            end
            # The electric energy term
            if epsilon != 0
                terms += epsilon, "Lsq", n + 1
                terms += 2 * l0 * epsilon, "L", n + 1
                terms += epsilon * l0^2, "Id", n + 1
            end
        end

        # The staggered mass term        
        if mu != 0
            terms += 0.5 * (-1)^nsite * mu, "Id", n
            terms += 0.5 * (-1)^nsite * mu, "Z", n
        end

        # The penalty term enforcing Gauss law
        if lambda != 0
            # We divide the constant part by N, so we can add it at every matter site, and in sum it is just 2λl_0^2
            terms += 2 * lambda * l0^2 / N, "Id", n
            # The Q^2 part of the penalty
            if isodd(nsite)
                terms += lambda, "Qoddsq", n
            else
                terms += lambda, "Qevensq", n
            end
            # L^2 on all the links
            if nsite < N
                terms += 2 * lambda, "Lsq", n + 1
            end
            # The terms -2L_n-1 L_n
            if nsite == 1
                terms -= 2 * lambda * l0, "L", n + 1
            elseif nsite == N
                terms -= 2 * lambda * l0, "L", n - 1
            else
                terms -= 2 * lambda, "L", n - 1, "L", n + 1
            end
            # The terms -2 Q_n Ln 
            if nsite == N
                if isodd(nsite)
                    terms -= 2 * lambda * l0, "Qodd", n
                else
                    terms -= 2 * lambda * l0, "Qeven", n
                end
            else
                if isodd(nsite)
                    terms -= 2 * lambda, "Qodd", n, "L", n + 1
                else
                    terms -= 2 * lambda, "Qeven", n, "L", n + 1
                end
            end
            # The terms 2 L_n-1 Q_n
            if nsite == 1
                terms += 2 * lambda * l0, "Qodd", n
            else
                if isodd(nsite)
                    terms += 2 * lambda, "L", n - 1, "Qodd", n
                else
                    terms += 2 * lambda, "L", n - 1, "Qeven", n
                end
            end
        end

        # The penalty term enforcing vanishing charge
        if eta != 0
            if isodd(nsite)
                terms += eta, "Qoddsq", n
            else
                terms += eta, "Qevensq", n
            end
            for k = n+2:2:length(s)
                nsiteprime = Int(round((k + 1) / 2))
                if isodd(nsite) && isodd(nsiteprime)
                    terms += 2 * eta, "Qodd", n, "Qodd", k
                elseif isodd(nsite) && iseven(nsiteprime)
                    terms += 2 * eta, "Qodd", n, "Qeven", k
                elseif iseven(nsite) && isodd(nsiteprime)
                    terms += 2 * eta, "Qeven", n, "Qodd", k
                else
                    terms += 2 * eta, "Qeven", n, "Qeven", k
                end
            end
        end
    end
    return MPO(terms, s)
end


############################################################################################
# Some observables
############################################################################################

"""
    getChargePenaltyMPO(s::Vector{Index{Int64}})

Get the MPO representation of the penalty term enforcing vanishing total charge
"""
function getChargePenaltyMPO(s::Vector{Index{Int64}})
    return getSchwingerHamiltonianMPO(0, 0, 0, 0, 0, 1, s)
end

"""
    getGaussLawPenaltyMPO(s::Vector{Index{Int64}})

Get the MPO representation of the penalty term enforcing vanishing Gauss law
"""
function getGaussLawPenaltyMPO(s::Vector{Index{Int64}})
    return getSchwingerHamiltonianMPO(0, 0, 0, 0, 1, 0, s)
end

"""
    getLMPO(n::Int, s::Vector{Index{Int64}})

Get the MPO for the electric field operator on the link between sites n and n+1
"""
function getLMPO(n::Int, s::Vector{Index{Int64}})
    terms = OpSum()
    n = 2*n
    terms += "L", n
    return MPO(terms, s)
end

############################################################################################
# MPS ground state calculations
############################################################################################

"""
    runSchwingerDMRG(N::Int, d::Int, J::Real, mu::Real, l0::Real, epsilon::Real, lambda::Real, eta::Real, weight::Real, nsweeps::Int=5000, mindim::Vector{Int64}=[20, 20, 50], maxdim::Vector{Int64}=[20, 40, 100], cutoff::Vector{<:Real}=[1E-11, 1E-11, 1E-11], noise::Vector{<:Real}=[1E-7, 1E-8, 0]; outputlevel::Int=1, etol::Real=1E-11, states::Vector{MPS}=Vector{MPS}(undef, 0), s::Vector{Index{Int64}}=Vector{Index{Int64}}(undef, 0))

Run an MPS computation for the Schwinger Hamiltonian running a two-site algorithm
"""
function runSchwingerDMRG(N::Int, d::Int, J::Real, mu::Real, l0::Real, epsilon::Real, lambda::Real, eta::Real, nsweeps::Int=100, mindim::Vector{Int64}=[5, 10, 10], maxdim::Vector{Int64}=[20, 60, 200], cutoff::Vector{<:Real}=[1E-11, 1E-11, 1E-12]; noise::Vector{<:Real}=[1E-7, 1E-8, 0])
    # Generate the sites
    s = getSites(N, d)
    # Generate the Hamiltonian MPO
    H = getSchwingerHamiltonianMPO(J, mu, epsilon, l0, lambda, eta, s)    
    #  Define the parameters for the different sweeps
    sweeps = Sweeps(nsweeps,
        [
            "maxdim" "mindim" "cutoff" "noise"
            maxdim mindim cutoff noise
        ])
    # Prepare a random wave function with D=5
    psi0 = randomMPS(s, 5)
    # Compute the ground state
    energy, psi = dmrg(H, psi0, sweeps, outputlevel=1)
    return energy, psi, s
end
