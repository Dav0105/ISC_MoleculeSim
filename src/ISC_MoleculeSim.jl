module ISC_MoleculeSim

using GLMakie, LinearAlgebra, Makie, CSV, TOML

export Molecule, Domain, generateSimulation

#region Structs
"""
Struct used to store Molecule properties  
`mass` = kg  
`radius` = m  
`chemicalFormula` = String   
`position` = (m, m)  
`speed` = m/s
"""
mutable struct Molecule
    # Mutable means that it can be modified
    chemicalFormula::String
    
    mass::Float64        # kg
    radius::Float64      # m

    position::Vector{Float64}
    speed::Vector{Float64}       # m/s

    pos_hist::Array{Vector{Float64}}
    speed_hist::Array{Vector{Float64}}
end

"""Represents the volume where the molecules will be."""
struct Domain
    # l_x::Float64 # m
    # l_y::Float64 # m
    # l_z::Float64 # m
    lims_x::Tuple{Float64, Float64} # m
    lims_y::Tuple{Float64, Float64} # m
    lims_z::Tuple{Float64, Float64} # m
end
#endregion

#region Computation functions
"""Returns the volume of the provided domain in m³."""
function getDomainVolume(d::Domain)::Number
    # return d.l_x * d.l_y * d.l_z
    return (d.lims_x[2] - d.lims_x[1]) * (d.lims_y[2] - d.lims_y[1]) * (d.lims_z[2] - d.lims_z[1])
end

"""Computes and applies the Movement of a single molecule with delta_time provided."""
function computeNextPosition(m::Molecule, dt::Number, g::Number)
    m.speed = m.speed .+ [0, 0, g * dt]
    m.position = m.position .+ (m.speed .* dt)
    return m.position
end

"""Function called when there is a collision. """
function computeCollisionVelocity(m1::Molecule, m2::Molecule)
    n::Vector = (m1.position .- m2.position) / sqrt(sum((m1.position .- m2.position).^2))
    # dot() = Produit Scalaire
    v1::Vector = m1.speed .- ((2 * m2.mass) / (m1.mass + m2.mass)) * dot((m1.speed .- m2.speed), n) .* n
    v2::Vector = m2.speed .+ ((2 * m1.mass) / (m1.mass + m2.mass)) * dot((m1.speed .- m2.speed), n) .* n

    m1.speed = v1
    m2.speed = v2
end

"""Checks if there are any collisions with provided molecules 
    and corrects the speed of molecules if there is collision."""
function computeMolsCollisions(mols::Array{Molecule})
    # Check collisions for this molecule with others
    for i in 1:length(mols)
        collides_with = nothing
        m = mols[i]
        for j in (i+1):length(mols)
            o_m = mols[j]

            if (detectMolsCollision(m, o_m))
                # Collision detected
                collides_with = o_m
                break
            end
        end

        if (collides_with !== nothing)
            # println("COLLISION DETECTED!!?!?!? $(m.chemicalFormula) => $(collides_with.chemicalFormula)")
            computeCollisionVelocity(m, collides_with)
        end
    end
end

"""Checks collisions between walls and provided molecules and corrects their
    velocity if needed"""
function checkWallCollisions(d::Domain, mols::Array{Molecule})::Nothing
    for m in mols
        mx, my, mz = m.position
        
        if (mx < d.lims_x[1]) || (mx > d.lims_x[2])
            m.speed[1] = -m.speed[1]
        end
        if (my < d.lims_y[1]) || (my > d.lims_y[2])
            m.speed[2] = -m.speed[2]
        end
        if (mz < d.lims_z[1]) || (mz > d.lims_z[2])
            m.speed[3] = -m.speed[3]
        end
    end
end

"""Compute all positions for all molecules until provided time, and returns an array containing [positions, speeds] for each t."""
function computePositions(mols::Array{Molecule}, domain::Domain, dt::Number, until::Number, g::Number, second_domain = nothing)

    for t in 0:dt:until
        # Changes domain at half time if second domain provided
        domain_to_use = (second_domain !== nothing && t >= until / 2) ? second_domain : domain

        # Check Collisions between molecules and walls
        checkWallCollisions(domain_to_use, mols)

        # Compute movement for each molecule
        for m in mols
            newPos = computeNextPosition(m, dt, g)
            push!(m.pos_hist,   newPos)      # current_m += newPos
            push!(m.speed_hist, m.speed)   # current_m += newSpeed
        end

        # Check Collisions between molecules
        computeMolsCollisions(mols)
    end

end

"""Returns True if m1 overlaps m2."""
function detectMolsCollision(m1::Molecule, m2::Molecule)::Bool
    if (m1 == m2) return false end

    r = m1.radius + m2.radius
    dist = sqrt(sum((m1.position .- m2.position).^2))

    return (dist <= r)
end
#endregion

#region Graphs generation

"""Displays a Makie graph containing the distributions of positions 
in Z of the molecules, at the end of the simulation."""
function generateGravityProbaGraph(
    mols::Array{Molecule}, 
    filename="./out/z_position_dist.png";
    fig_size::Tuple{Int, Int} = (1000, 800)
)
    last_pos_hist = []
    for mol in mols
        append!(last_pos_hist, last(mol.pos_hist)[3])
    end

    f = Figure(size = fig_size)
    ax = Axis(f[1, 1],
        title = "Distribution of molecules positions in Z at the end of the simulation",
        xlabel = "Position in Z",
        ylabel = "Number of molecules",
    )
    hist!(ax, last_pos_hist)
    save(filename, f)           # Save figure
    # display(f)
    return
end

"""Displays a Makie graph of the average m*v^2 across time for the provided molecules."""
function generateMeanMv2Graph(
    mols::Array{Molecule},
    delta_t::Number,
    filename::String = "./out/mv2_vs_time.png";
    fig_size::Tuple{Int, Int} = (1000, 800)
)
    if isempty(mols) || isempty(first(mols).speed_hist)
        error("Cannot plot <mv^2> without simulated molecule speed history")
    end

    frame_count = length(first(mols).speed_hist)
    times = collect(1:frame_count) .* delta_t
    mean_mv2 = zeros(frame_count)

    for frame in 1:frame_count
        total_mv2 = 0.0
        for mol in mols
            speed = mol.speed_hist[frame]
            total_mv2 += mol.mass * dot(speed, speed)
        end
        mean_mv2[frame] = total_mv2 / length(mols)
    end

    f = Figure(size = fig_size)
    ax = Axis(
        f[1, 1],
        title = "Average <m v^2> across time",
        xlabel = "Time (s)",
        ylabel = "<m v^2> (kg·m²/s²)",
    )
    lines!(ax, times, mean_mv2)
    save(filename, f)
    # display(f)
    return
end

"""Displays a Makie graph of the average temperature across time for the provided molecules."""
function generateTemperatureGraph(
    mols::Array{Molecule},
    delta_t::Number,
    filename::String = "./out/temperature_vs_time.png";
    fig_size::Tuple{Int, Int} = (1000, 800)
)
    boltzmann_constant = 1.380649 * 10^-23 # J/K

    frame_count = length(first(mols).speed_hist)
    mean_mv2 = zeros(frame_count)
    times = collect(1:frame_count) .* delta_t
    for frame in 1:frame_count

        # Compute m * v^2 each molecule at this frame
        for mol in mols
            speed = mol.speed_hist[frame]
            mean_mv2[frame] += mol.mass * norm(speed)^2
        end

        # Average of m * v^2 over number of molecules
        mean_mv2[frame] /= length(mols)

        # Convert to temperature using T = (m * <v^2>) / (3 * k_B)
        mean_mv2[frame] /= (3 * boltzmann_constant)
    end

    # PLOT THAT
    f = Figure(size = fig_size)
    ax = Axis(
        f[1, 1],
        title = "Temperature (proportional to m·<v²>) across time for each molecule type",
        xlabel = "Time (s)",
        ylabel = "Temperature (K)",
    )
    lines!(ax, times, mean_mv2)

    save(filename, f)
    # display(f)
    return
end

"""Displays a Makie graph of the probability density of molecules positions in Z."""
function generateZProbabilityGraph(
    mols::Array{Molecule},
    num_bins::Int = 10,
    filename::String = "./out/z_position_probability.png";
    fig_size::Tuple{Int, Int} = (1000, 800)
)
    z_positions = [last(mol.pos_hist)[3] for mol in mols]

    f = Figure(size = fig_size)
    ax = Axis(f[1, 1],
        title = "Distribution of molecules positions in Z at the end of the simulation",
        xlabel = "Position in Z",
        ylabel = "Probability",
    )
    hist!(ax, z_positions, bins=num_bins, normalization=:probability)
    save(filename, f)  # Save figure
    # display(f)
    return
end

function generatePressureGraph(
    mols::Array{Molecule},
    delta_t::Number,
    domain::Domain,
    filename::String = "./out/pressure_vs_time.png";
    fig_size::Tuple{Int, Int} = (1000, 800)
)

    frame_count = length(first(mols).speed_hist)
    mean_mv2 = zeros(frame_count)
    times = collect(1:frame_count) .* delta_t
    for frame in 1:frame_count

        # Compute m * v^2 each molecule at this frame
        for mol in mols
            speed = mol.speed_hist[frame]
            mean_mv2[frame] += mol.mass * norm(speed)^2
        end

        # Average of m * v^2 over number of molecules
        mean_mv2[frame] /= length(mols)

        # Convert to Pressure using P = (1/3) * (N/V) * <mv^2>
        mean_mv2[frame] /= (3 * getDomainVolume(domain))
    end

    # PLOT THAT
    f = Figure(size = fig_size)
    ax = Axis(
        f[1, 1],
        title = "Pressure (proportional to m·<v²>) across time for each molecule type",
        xlabel = "Time (s)",
        ylabel = "Pressure (Pa)",
    )
    lines!(ax, times, mean_mv2)

    save(filename, f)
    # display(f)
    return
end

"""Displays a Makie graph of the entropy of the system across time."""
function generateEntropyGraph(
    mols::Array{Molecule},
    domain::Domain,
    delta_t::Number,
    filename::String = "./out/entropy_vs_time.png";
    fig_size::Tuple{Int, Int} = (1000, 800),
    probability_bins::Tuple{Int, Int, Int, Int} = (10, 10, 10, 200),
    second_domain = nothing
)
    # Get entropy for each frame
    if (second_domain !== nothing)
        entropies = calculateEntropies(mols, second_domain, probability_bins)
    else
        entropies = calculateEntropies(mols, domain, probability_bins)
    end

    # To display time in x axis
    frame_count = length(first(mols).speed_hist)
    times = collect(1:frame_count) .* delta_t

    # PLOT THIS :D
    f = Figure(size = fig_size)
    ax = Axis(
        f[1, 1],
        title = "Entropy across time for the system of molecules",
        xlabel = "Time (s)",
        ylabel = "Entropy (bits)",
    )

    lines!(ax, times, entropies)
    save(filename, f)
    # display(f)
    return
end
#endregion

#region Entropy calculation
"""Returns an array containing the probability distribution of positions and velocities of the provided molecules for each time."""
function getProbabilityDistribution(mols::Array{Molecule}, domain::Domain, probability_bins::Tuple{Int, Int, Int, Int} = (10, 10, 10, 200))

    function getProba(elements::Array{Float64}, bin_domain::Tuple{Float64, Float64}, num_bins::Int = 10)
        z_min, z_max = bin_domain
        bin_width = (z_max - z_min) / num_bins

        # Count the number of positions in each bin
        counts = zeros(Int, num_bins)
        for e in elements
            for bin_idx in 1:num_bins
                curr_min_z = z_min + (bin_idx - 1) * bin_width
                curr_max_z = z_min + bin_idx * bin_width
                if curr_min_z <= e < curr_max_z
                    counts[bin_idx] += 1
                    break
                end
            end
        end
        
        # Normalize to get probabilities
        probabilities = [c / length(elements) for c in counts]

        return probabilities
    end

    x_probas = []; y_probas = []; z_probas = []; v2_probas = []
    for t in 1:length(first(mols).pos_hist)        
        push!(x_probas, getProba([mol.pos_hist[t][1] for mol in mols], (domain.lims_x[1], domain.lims_x[2]), probability_bins[1]))
        push!(y_probas, getProba([mol.pos_hist[t][2] for mol in mols], (domain.lims_y[1], domain.lims_y[2]), probability_bins[2]))
        push!(z_probas, getProba([mol.pos_hist[t][3] for mol in mols], (domain.lims_z[1], domain.lims_z[2]), probability_bins[3]))
        push!(v2_probas, getProba([norm(mol.speed_hist[t]) for mol in mols], (0.0, 200_000.0), probability_bins[4]))
    end

    return x_probas, y_probas, z_probas, v2_probas
end

"""Returns an array containing the entropy of the system at each time.  
The entropy is calculated using the formula `H(X, Y, Z, < v^2 >) = H(X) + H(Y) + H(Z) + H(< v^2 >)`."""
function calculateEntropies(mols::Array{Molecule}, domain::Domain, probability_bins::Tuple{Int, Int, Int, Int} = (10, 10, 10, 200))
    # Compute the probability distribution of positions and velocities
    pxs, pys, pzs, pv2s = getProbabilityDistribution(mols, domain, probability_bins)

    # H(X) = -Σ P(x) log P(x)
    function entropy(probas)
        # iszero() to avoid NaN values when proba is 0.0
        return -sum(p -> iszero(p) ? 0.0 : p * log2(p), probas)
    end

    # Get entropy for each probability distribution
    # H(X, Y, Z, < v^2 >) = H(X) + H(Y) + H(Z) + H(< v^2 >)
    entropies = [entropy(px) + entropy(py) + entropy(pz) + entropy(pv2) for (px, py, pz, pv2) in zip(pxs, pys, pzs, pv2s)]

    return entropies
end
#endregion

#region Main function
"""Generate simulation with provided settings and outputs it to `./out` folder.

**Params:**  
domain : Domain of the simulation  
mols : Array of Molecules to add to the simulation  
delta_t : Time step used for the simulation (in seconds)  
until : Time until which the simulation will be generated (in seconds)  
framerate : Framerate of the output video (in frames per second)  
  
framestep : Number of frames to skip between each frame of the output video  
export_to_csv : Whether to export the simulation data to a CSV file (default: true) (not implemented yet, only exports settings to a TOML file for now)  
output_path : Path where the output video and graphs will be saved (without extension, default: "./out/animation")  
g : Gravitational acceleration to apply to molecules (in m/s², default: -9.81)  
second_domain : Optional second domain to use for the second half of the simulation (default: nothing)  
probability_bins : Number of bins to use for the probability distribution when calculating entropy (x, y, z and v^2)
"""
function generateSimulation(
    domain::Domain, mols::Array{Molecule}, delta_t::Number, until::Number, framerate::Int; 
    framestep::Int= 30, exportToCSV::Bool = false, output_path::String = "./out/animation", g::Number=-9.81,
    second_domain = nothing, probability_bins::Tuple{Int, Int, Int, Int} = (10, 10, 10, 200)
)
    # Domain volume
    println("Domain volume : " * string(getDomainVolume(domain)) * " m³")
    if (second_domain !== nothing)
        println("Second domain volume : " * string(getDomainVolume(second_domain)) * " m³")
    end

    timestamps = 0:delta_t:until
    total_frames = length(timestamps)

    # Observable indicating current frame being generated
    # (Used in the record section)
    frame = Observable(1)
    # Observable containing time info
    T = Observable(0.0)

    @time "Time to generate positions" begin
        # Compute all positions for all molecules
        computePositions(mols, domain, delta_t, until, g, second_domain)
    end

    # Separate molecules into different arrays by their chemical formula (for color purposes)
    """mols_by_formula = {
        "He" => [m1, m2, m3],
        "Ar" => [m4, m5],
        ...
    }"""
    mols_by_formula = Dict{String, Array{Molecule}}()
    for m in mols
        if haskey(mols_by_formula, m.chemicalFormula)
            push!(mols_by_formula[m.chemicalFormula], m)
        else
            mols_by_formula[m.chemicalFormula] = [m]
        end
    end

    # Prepare figure
    fig = Figure(size = (800, 600))

    if second_domain !== nothing
        xlims = (second_domain.lims_x[1], second_domain.lims_x[2])
        ylims = (second_domain.lims_y[1], second_domain.lims_y[2])
        zlims = (second_domain.lims_z[1], second_domain.lims_z[2])
    else
        xlims = (domain.lims_x[1], domain.lims_x[2])
        ylims = (domain.lims_y[1], domain.lims_y[2])
        zlims = (domain.lims_z[1], domain.lims_z[2])
    end

    aspect_x = xlims[2] - xlims[1]
    aspect_y = ylims[2] - ylims[1]
    aspect_z = zlims[2] - zlims[1]

    ax = Axis3(
        fig[1, 1], 
        # perspectiveness = 0.5,
        aspect = (aspect_x, aspect_y, aspect_z), 
        title = @lift("t = $($T) s"),
        limits=(xlims, ylims, zlims),
        azimuth = 0.3 * pi
    )    
    
    # Init 1st position for scatters (1 scatter per molecule type)
    scatter_plots = []
    for i in 1:length(mols_by_formula)
        xs, ys, zs = mols[1].pos_hist[1]
        # println(mol_formula)
        # println(length(mols_list))
        sc_plot = scatter!(ax, xs, ys, zs, markersize=10, color=i, colorrange = (1, length(mols_by_formula)))
        push!(scatter_plots, sc_plot)
    end
    
    @time "Time to generate video" record(fig, output_path * ".mp4", 1:framestep:total_frames; framerate = framerate) do f
        frame[] = f
        T[] = timestamps[f]

        # Generate a scatter for each mol type (each having their color)
        for (i, mol_formula) in enumerate(keys(mols_by_formula))
            # Update scatter with current positions
            xs = []; ys = []; zs = []
            
            mols_curr_form = mols_by_formula[mol_formula]
            for m in mols_curr_form     # foreach molecule
                x, y, z = m.pos_hist[f] # Get positions at frame f for molecule m
                push!(xs, x)
                push!(ys, y)
                push!(zs, z)
            end

            scatter_plots[i][1] = xs  # xs = [mol1_x, mol2_x, mol3_z] at frame f
            scatter_plots[i][2] = ys  # y
            scatter_plots[i][3] = zs  # z
        end

    end

    # Exports simultation settings in a TOML file
    # And CSV file containing pos and speeds (not done yet)
    if (exportToCSV)
        # CSV Structure :
        # frame, speed_x, speed_y, speed_z, pos_x, pos_y, pos_z
        # 1, [dataForMol1]
        # 1, [dataForMol2]
        # ...

        # CSV.write(output_path * "_speed.csv", df)

        # Export settings
        settings = Dict( 
            "simulation" => Dict(
                "delta_t" => delta_t,
                "until" => until,
                "domain_x" => [domain.lims_x[1], domain.lims_x[2]],
                "domain_y" => [domain.lims_y[1], domain.lims_y[2]],
                "domain_z" => [domain.lims_z[1], domain.lims_z[2]],
                "framerate" => framerate,
                "g" => g
            )
        )
        open(output_path * ".toml", "w") do io
            TOML.print(io, settings)
        end
    end


    println("Video saved in " * output_path * " ! :)")

    # Graph generation
    @time "Time to generate graphs" begin
        size = (1000, 800)
        generateGravityProbaGraph(mols, output_path * "_z_positions_dist.png", fig_size=size)
        generateMeanMv2Graph(mols, delta_t, output_path * "_mv2_vs_time.png", fig_size=size)
        generateTemperatureGraph(mols, delta_t, output_path * "_temperature_vs_time.png", fig_size=size)
        generatePressureGraph(mols, delta_t, domain,output_path * "_pressure_vs_time.png", fig_size=size)
        generateZProbabilityGraph(mols, 10, output_path * "_z_position_probability.png", fig_size=size)
        generateEntropyGraph(mols, domain, delta_t, output_path * "_entropy_vs_time.png", fig_size=size, probability_bins=probability_bins, second_domain=second_domain)
    end

end
#endregion

end # module ISC_MoleculeSim
