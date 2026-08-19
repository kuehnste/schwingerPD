##############################################################
# Driver to run MPS calculations
##############################################################
using DelimitedFiles
using ITensors, ITensorMPS

include("../src/SchwingerModelGaugeFieldIntegrated.jl")


let
    ##################################################################
    # Parameters of the model
    ##################################################################

    # Number of Spins
    N = 20
    # The number of steps we take scan the background field
    nsteps = 20

    # Inverse lattice spacing squared in units of the coupling
    x = 1.0
    # Mass in units of the coupling
    mg = 0.5
    # Maximum background electric field
    l0max = 2.0
    # Strength of the penalty term
    eta = 100.0

    # Convert x and m/g in the corresponding parameters needed for the Schwinger Hamiltonian
    J = x
    mu = 2 * mg * sqrt(x)
    epsilon = 1.0
    # Steps in the background field
    dl0 = l0max/(nsteps-1)
    # Run the calculation
    count = 1
    results = zeros(nsteps, N+1)
    for l0 in 0:dl0:l0max
        println("--> Working on l0 = ", l0)
        energy, gs, s = runSchwingerDMRG(N, J, mu, l0, epsilon, eta)
        # Compute the expected value of Pauli-Z everywhere
        Zvals = expect(gs, "Z")
        # Get the electric field
        Lvals = getElectricField(Zvals)
        # Compute the expected value for the total charge and check if we have vanishing total charge
        ChargePenaltyMPO = getChargePenaltyMPO(s)
        P = inner(gs', ChargePenaltyMPO, gs)
        if abs(P)>1E-10
            @warn string("Penalty enforcing vanishing total charge has a value of ", P)
        end
        results[count, :] .= [l0; energy; Lvals]
        count += 1
    end
    # Save the results to a file
    fname = string("schwinger_gauge_field_integrated_N", N, "_x", x, "_mg", mg, "_eta", eta, ".txt")
    open(fname, "w") do io
        writedlm(io, results)
    end
end
