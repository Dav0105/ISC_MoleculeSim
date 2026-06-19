module ISC_LennardJonesMolSim

using Distributions, GLMakie, LinearAlgebra

#region Structs n consts
"""Struct used to store Molecule properties
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

    force_to_apply::Vector{Float64} # N

    pos_hist::Array{Vector{Float64}}
    speed_hist::Array{Vector{Float64}}
end

"""Represents the volume where the molecules will be.
lims_ = (min, max) for A in {x, y}"""
struct Domain
    # l_x::Float64 # m
    # l_y::Float64 # m
    lims_x::Tuple{Float64,Float64} # m
    lims_y::Tuple{Float64,Float64} # m
end

const boltzmann_constant = 1.380649 * 10^-23 # J/K
const u = 1.660_538_921 * (10^-27) # [kg]

#endregion

#region Helpers
"""Calculate the Lennard-Jones potential between two molecules"
`r` = distance between the two molecules (m)
`epsilon` (ε) = depth of the potential well (J)
`sigma` (σ) = finite distance at which the inter-particle potential is zero (m)
"""
function lennard_jones_potential(r::Float64, epsilon::Float64, sigma::Float64)
    return 4 * epsilon * ((sigma / r)^12 - (sigma / r)^6)
end

"""Calculate the Lennard-Jones force between two molecules
`r_i` = position of molecule i (m)
`r_j` = position of molecule j (m)
`epsilon` (ε) = depth of the potential well (J)
`sigma` (σ) = finite distance at which the inter-particle potential is zero (m)
"""
function lennard_jones_force(
    r_i::Vector{Float64}, r_j::Vector{Float64},
    epsilon::Float64, sigma::Float64
)::Vector{Float64}
    r_ij = r_j - r_i
    r = norm(r_ij)

    # Prevent division by zero and singularities, but keep the short-range
    # repulsion so overlapping particles are pushed apart instead of merging.
    if !isfinite(r) || r >= 3 * sigma
        return zeros(2)
    end

    # r = max(r, 0.5 * sigma)

    return 24 * epsilon * (1 / r^2) * ((2 * (sigma / r)^12) - (sigma / r)^6) * r_ij
end

function calculate_force_for_molecule(mol::Molecule, molecules::Vector{Molecule}, epsilon::Float64, sigma::Float64)
    total_force::Vector{Float64} = zeros(2)  # Initialize total force as a 2D vector

    for other_mol in molecules
        if other_mol !== mol  # Avoid calculating force with itself
            force = lennard_jones_force(other_mol.position, mol.position, epsilon, sigma)
            total_force .+= force
        end
    end

    return total_force

end

function calculate_system_temp(mols::Vector{Molecule})::Float64
    total_mv2::Float64 = 0.0
    # Calculate avg mv^2 for each frame
    for mol in mols
        total_mv2 += mol.mass * norm(mol.speed)^2
    end

    avg_mv2 = total_mv2 / length(mols)
    res = avg_mv2 / (2 * boltzmann_constant)

    # Guard against NaN
    if isnan(res) || res < 0
        print("Warning: Invalid temperature calculation. Returning 0.0 K\n")
        return 0.0
    end

    # print("Curr temp: ", res, " K\n")
    return res
end

function rescale_velocities(mols::Vector{Molecule}, target_temp::Float64)::Float64
    current_temp = calculate_system_temp(mols)

    rescale_factor = get_rescale_velocity_factor(current_temp, target_temp)

    for mol in mols
        mol.speed *= rescale_factor
    end

    return current_temp
end

function get_rescale_velocity_factor(current_temp::Float64, target_temp::Float64)
    scaling_factor = sqrt(target_temp / current_temp)
    return scaling_factor
end

function computeNextPosition(mol::Molecule, force::Vector{Float64}, dt::Float64)
    acceleration = force ./ mol.mass

    # Euler integration
    mol.speed = mol.speed .+ (acceleration * dt)
    mol.position = mol.position .+ (mol.speed .* dt)

    # Store history
    push!(mol.pos_hist, copy(mol.position))
    push!(mol.speed_hist, copy(mol.speed))
end

"""Checks if the molecule is within the domain, if not, it reflects the speed in the direction of the boundary"""
function checkDomainBounds(mol::Molecule, domain::Domain)
    mx, my = mol.position
    
    # Check x bounds
    if (mx < domain.lims_x[1]) || (mx > domain.lims_x[2])
        if mx < domain.lims_x[1]
            dist_domain_mol = domain.lims_x[1] - mx
            mol.position[1] = mol.position[1] + 2 * dist_domain_mol
        elseif mx > domain.lims_x[2]
            dist_domain_mol = domain.lims_x[2] - mx
            mol.position[1] = mol.position[1] - 2 * dist_domain_mol
        end
        mol.speed[1] = -mol.speed[1]
    end

    # Check y bounds
    if (my < domain.lims_y[1]) || (my > domain.lims_y[2])
        if my < domain.lims_y[1]
            dist_domain_mol = domain.lims_y[1] - my
            mol.position[2] = mol.position[2] + 2 * dist_domain_mol
        elseif my > domain.lims_y[2]
            dist_domain_mol = domain.lims_y[2] - my
            mol.position[2] = mol.position[2] - 2 * dist_domain_mol
        end
        mol.speed[2] = -mol.speed[2]
    end

end


function computePositions(mols::Vector{Molecule}, domain::Domain, dt::Float64, until::Float64, epsilon::Float64, sigma::Float64, target_temp::Float64=10.0)::Vector{Float64}
    temperatures = []
    
    for t in 0:dt:until
        # Calculate forces for mols
        for mol in mols
            mol.force_to_apply = calculate_force_for_molecule(mol, mols, epsilon, sigma)
        end
        # Temporal integration (Euler method)
        for mol in mols
            computeNextPosition(mol, mol.force_to_apply, dt)

            # Specular reflection
            checkDomainBounds(mol, domain)
        end
        # Velocity rescaling
        current_temp::Float64 = rescale_velocities(mols, target_temp)
        push!(temperatures, current_temp)
    end

    return temperatures
end

function generateTempGraph(temps::Vector{Float64}, delta_t::Float64, until::Float64, output_path::String="./out/temperature_graph.png")
    timestamps = 0:delta_t:until
    fig = Figure()
    ax = Axis(
        fig[1, 1], 
        title="Temperature vs Time", 
        xlabel="Time (s)", 
        ylabel="Temperature (K)"
    )
    lines!(ax, timestamps, temps, color=:blue)
    save(output_path, fig)
end

#endregion

function generateSimulation(
    domain::Domain, mols::Vector{Molecule}, delta_t::Number, until::Number, framerate::Int, epsilon::Float64=1.0, sigma::Float64=1.0, target_temp=10.0;
    framestep::Int=30, output_path::String="./out/2D_animation"
)
    timestamps = 0:delta_t:until
    total_frames = length(timestamps)

    # Observable indicating current frame being generated
    # (Used in the record section)
    frame = Observable(1)
    # Observable containing time info
    T = Observable(0.0)

    curr_temp = Observable(0.0)

    # Initial calculated temperature
    initial_temp = calculate_system_temp(mols)
    println("Initial temperature: $initial_temp K")

    @time "Time to generate positions" begin
        # Compute all positions for all molecules
        temperatures = computePositions(mols, domain, delta_t, until, epsilon, sigma, target_temp)
    end

    xlims = (domain.lims_x[1], domain.lims_x[2])
    ylims = (domain.lims_y[1], domain.lims_y[2])

    fig = Figure(size=(800, 800))
    ax = Axis(
        fig[1, 1],
        # perspectiveness = 0.5,
        # aspect=(aspect_x, aspect_y),
        title=@lift("t = $($T) s"),
        subtitle=@lift("Tcurr = $($curr_temp) K"),
        limits=(xlims, ylims),
    )

    # Init 1st positions
    xs, ys = mols[1].pos_hist[1]
    sc_plot = scatter!(ax, xs, ys, markersize=13)


    @time "Time to generate video" record(fig, output_path * ".mp4", 1:framestep:total_frames; framerate=framerate) do f
        frame[] = f
        T[] = timestamps[f]
        curr_temp[] = temperatures[f]

        xs = Float64[]
        ys = Float64[]

        for (i, mol) in enumerate(mols)
            x, y = mol.pos_hist[f] # Get positions at frame f for molecule m
            push!(xs, x)
            push!(ys, y)
        end
        
        sc_plot[1] = xs  # xs = [mol1_x, mol2_x, mol3_z] at frame f
        sc_plot[2] = ys  # y
    end

    println("Video saved in " * output_path * " ! :)")

    generateTempGraph(temperatures, delta_t, until, "./out/temperature_graph.png")
end

#region Main
function generate_random_positions(domain::Domain, existing_mols::Vector{Molecule}, min_distance::Float64)
    while true
        pos = [
            rand() * (domain.lims_x[2] - domain.lims_x[1]) + domain.lims_x[1],
            rand() * (domain.lims_y[2] - domain.lims_y[1]) + domain.lims_y[1]
        ]
        
        # Check distance from all existing molecules
        too_close = false
        for mol in existing_mols
            if norm(pos - mol.position) < min_distance
                too_close = true
                break
            end
        end
        
        if !too_close
            return pos
        end
    end
end

function main()
    domain::Domain = Domain(
        (-5 * 10^-9, 5 * 10^-9),
        (-5 * 10^-9, 5 * 10^-9)
    )

    delta_t = 1 * 10^-15 # s
    until = 5 * 10^-11 # s

    num_mols = 100
    T = 40.0 # K
    sigma = 2.74 * 10^-10 # m
    epsilon = 4.91511044 * 10^-22 # J


    mols::Vector{Molecule} = []
    mol_mass = 20.1797 * u
    for i in 1:num_mols
        position = generate_random_positions(domain, mols, 0.9 * sigma)
        speed = rand(Normal(0, sqrt(boltzmann_constant * T / mol_mass)), 2)
        push!(mols, Molecule(
            "Ne", mol_mass, 38 * 10^-12, position, speed, zeros(2), [], []
        ))
    end

    generateSimulation(domain, mols, delta_t, until, 30, epsilon, sigma, T; framestep=100)
end
#endregion

main()

end # module ISC_LennardJonesMolSim
