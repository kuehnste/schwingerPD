##############################################################
# Visualize the results from the MPS calculations
##############################################################
using DelimitedFiles
using LaTeXStrings
using PyPlot


"""
    plotEnergyElectricField(fname::Vector{String}, naverage::Int=2)

Plot the energy and the average electric field in the center of the system as a function of the background field. 
The input naverage indicates over how many fields the electric field will be averaged
"""
function plotEnergyElectricField(fnames::Vector{String}, naverage::Int=2, labelstring::Vector{<:AbstractString}=Vector{String}(undef, 0))
    # Prepare a figure
    h = figure(figsize=(17, 6))
    rc("text", usetex=true)
    rc("font", family="serif", weight="bold", size="20")
    rc("axes", linewidth=2)
    rc("xtick.major", size=10)
    rc("xtick.major", width=2)
    rc("xtick.minor", size=5)
    rc("xtick.minor", width=2)
    rc("ytick.major", size=10)
    rc("ytick.major", width=2)
    rc("ytick.minor", size=5)
    rc("ytick.minor", width=2)

    # The markers we use
    markers = ["s"; "^"; "o"; "h"; "<"; "p"; ">"; "v"; "+"; "H"; "x"; "+"; "."; "P"]
    # Markersize
    msize = 10
    # Line width
    lwidth = 3

    # Prepare two panels, one for the energy, one for the electric field
    ax1 = plt.subplot2grid((1, 2), (0, 0))
    ax2 = plt.subplot2grid((1, 2), (0, 1))

    # Read the data and plot
    count = 1
    for file in fnames
        if !isfile(file)
            @warn string("File ", file, " not found")
        else
            lstring = file
            data = readdlm(file)
            # Extract the system size
            N = size(data, 2) - 1
            # Average the electric field in the center
            ind = Int(round(N/2))
            Laverage = sum(data[:, (2+ind):(2+ind+naverage-1)], dims=2)
            # In case custom label strings are provided, we use these
            if !isempty(labelstring)
                lstring = labelstring[count]
            end
            # Plot the quantities
            ax1.plot(data[:, 1], data[:, 2], ":"*markers[count], ms=msize, lw=lwidth, label=lstring)
            ax2.plot(data[:, 1], Laverage, ":"*markers[count], ms=msize, lw=lwidth, label=lstring)
            if count<length(markers)
                count += 1
            end
        end
    end
    ax1.set_xlabel(L"l_0 = \theta/2\pi")
    ax1.set_ylabel(L"E_0")
    ax1.set_title("Ground state energy")
    ax1.tick_params(which="both", direction="in")

    ax2.set_xlabel(L"l_0 = \theta/2\pi")
    ax2.set_ylabel(L"\sum_{n\in N/2}\langle L_n\rangle")
    ax2.set_title("Average electric field")
    ax2.legend()
    ax2.tick_params(which="both", direction="in")
end