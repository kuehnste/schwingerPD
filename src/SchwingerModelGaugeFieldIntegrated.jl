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
function getSites(N::Int)
    return siteinds("S=1/2", N)
end


############################################################################################
# The Hamiltonian
############################################################################################

"""
    getSchwingerHamiltonianMPO(J::Real, mu::Real, epsilon::Real, l0::Real, eta::Real, eta::Real, s::Vector{Index{Int64}}, q::Vector{Int}=Vector{Int64}(undef, 0))

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
function getSchwingerHamiltonianMPO(J::Real, mu::Real, epsilon::Real, l0::Real, eta::Real, s::Vector{Index{Int64}})
    # Build the Hamiltonian MPO, seemingly unnecessary if-statements, are there to prevent terms in the MPO with prefactors that are numerical zeros and could potentially increase the bond dimension without having any effect            
    terms = OpSum()

    N = length(s)

    # Precompute the site-dependent numerical factors arising form the electric energy term
    prefactor = zeros(N-1)
    prefactor[N-1] = -1
    for k = (N-2):-1:1
        prefactor[k] = prefactor[k+1]
        if isodd(k)
            prefactor[k]-=1
        end
    end

    prefactor2=(0.25-l0)*N/2
    for n = 1:(N-1)
        prefactor2 += (0.25 * (N-n))
    end

    terms += prefactor2*epsilon, "Id", 1
    # Constant term from the elctric field energy
    terms += (N-1)*l0^2*epsilon, "Id", 1

    for n = 1:length(s)
        #terms += 1/2, "Id",n # constant term from expanding L^2
        if n < length(s)
            # The kinetic term
            if J != 0
                terms -= J, "S+", n, "S-", n + 1
                terms -= J, "S-", n, "S+", n + 1
            end
            # The electric energy term
            terms += l0*(N-n)*epsilon + 0.5*prefactor[n]*epsilon, "Z", n            
            for m = 1:(n-1)
                terms += 0.5 * (N-n) * epsilon, "Z", n, "Z", m
            end
        end

        # The staggered mass term        
        if mu != 0
            terms += 0.5 * (-1)^n * mu, "Id", n
            terms += 0.5 * (-1)^n * mu, "Z", n
        end

        # The penalty term enforcing vanishing charge
        if eta != 0
            if isodd(n)
                terms += eta, "Qoddsq", n
            else
                terms += eta, "Qevensq", n
            end
            for k = (n+1):length(s)
                if isodd(n) && isodd(k)
                    terms += 2 * eta, "Qodd", n, "Qodd", k
                elseif isodd(n) && iseven(k)
                    terms += 2 * eta, "Qodd", n, "Qeven", k
                elseif iseven(n) && isodd(k)
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
# Computing observables
############################################################################################

"""
    getChargePenaltyMPO(s::Vector{Index{Int64}})

Get the MPO representation of the penalty term enforcing vanishing total charge
"""
function getChargePenaltyMPO(s::Vector{Index{Int64}})
    return getSchwingerHamiltonianMPO(0, 0, 0, 0, 1, s)
end

"""
    getElectricField(Zvalues::Vector{<:Number})

Given a vector containing the expected value of Pauli-Z on each site, reconstruct the corresponding electric field on the links
"""
function getElectricField(Zvalues::Vector{<:Number})
    # Extract the system size
    N = length(Zvalues)
    # From the spin expectation value reconstruct the charge
    Qvalues = 0.5*Zvalues .+ (-1) .^ collect(1:N)
    # The electric field can be obtained as cumulative sum of the charge values, where we consider there is no nontrivial link to the right of the last matter site
    Lvalues = cumsum(Qvalues[1:(end-1)])
    return Lvalues
end


############################################################################################
# MPS ground state calculations
############################################################################################

"""
    runSchwingerDMRG(N::Int, J::Real, mu::Real, l0::Real, epsilon::Real, eta::Real, nsweeps::Int=100, mindim::Vector{Int64}=[5, 10, 10], maxdim::Vector{Int64}=[20, 60, 200], cutoff::Vector{<:Real}=[1E-11, 1E-11, 1E-12]; noise::Vector{<:Real}=[1E-7, 1E-8, 0])

Run an MPS computation for the Schwinger Hamiltonian running a two-site MPS algorithm for variational ground state search
"""
function runSchwingerDMRG(N::Int, J::Real, mu::Real, l0::Real, epsilon::Real, eta::Real, nsweeps::Int=100, mindim::Vector{Int64}=[5, 10, 10], maxdim::Vector{Int64}=[20, 60, 200], cutoff::Vector{<:Real}=[1E-11, 1E-11, 1E-12]; noise::Vector{<:Real}=[1E-7, 1E-8, 0])
    # Generate the sites
    s = getSites(N)
    # Generate the Hamiltonian MPO
    H = getSchwingerHamiltonianMPO(J, mu, epsilon, l0, eta, s)
    # Define the parameters for the different sweeps
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
