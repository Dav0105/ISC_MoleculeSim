module ISC_MoleculeSim

using GLMakie, LinearAlgebra, Makie, CSV, TOML

export Molecule, Domain, generateSimulation

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
    
    mass::Number        # kg
    radius::Number      # m

    position::Vector
    speed::Vector       # m/s
end

"""Represents the volume where the molecules will be."""
struct Domain
    l_x::Number # m
    l_y::Number # m
    l_z::Number # m
end

"""Returns the volume of the provided domain in m³."""
function getDomainVolume(d::Domain)::Number
    return d.l_x * d.l_y * d.l_z
end

"""Computes the Movement of a molecule with delta_time provided."""
function computeMovement(m::Molecule, dt::Number)
    m.position = m.position .+ (m.speed * dt)
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
        dlx = d.l_x / 2
        dly = d.l_y / 2
        dlz = d.l_z / 2
        
        if (mx < -dlx) || (mx > dlx)
            m.speed[1] = -m.speed[1]
        end
        if (my < -dly) || (my > dly)
            m.speed[2] = -m.speed[2]
        end
        if (mz < -dlz) || (mz > dlz)
            m.speed[3] = -m.speed[3]
        end
    end
end

"""Compute all positions for all molecules until provided time, and returns an array containing [positions, speeds] for each t."""
function computePositions(mols::Array{Molecule}, domain::Domain, dt::Number, until::Number)::Tuple{Vector, Vector}
    pos::Vector{Vector} = []
    speed::Vector{Vector} = []

    # Create new array for molecule if doesn't exist
    for _ in 1:length(mols)
        push!(pos, [])
        push!(speed, [])
    end

    for t in 0:dt:until
        # Check Collisions between molecules and walls
        checkWallCollisions(domain, mols)

        # Compute movement for each molecule
        for (idx, m) in enumerate(mols)
            newPos = computeMovement(m, dt)
            push!(pos[idx], newPos)      # current_m += newPos
            push!(speed[idx], m.speed)   # current_m += newSpeed
        end

        # Check Collisions between molecules
        computeMolsCollisions(mols)


    end

    return (pos, speed)
end

"""Returns True if m1 overlaps m2."""
function detectMolsCollision(m1::Molecule, m2::Molecule)::Bool
    if (m1 == m2) return false end

    r = m1.radius + m2.radius
    dist = sqrt(sum((m1.position .- m2.position).^2))

    return (dist <= r)
end

"""Generate simulation with provided settings and outputs it to `./out` folder."""
function generateSimulation(domain::Domain, mols::Array{Molecule}, delta_t::Number, until::Number, framerate::Int, exportToCSV::Bool = true)::Nothing
    output_path = "out/animation.mp4"
    
    println("Domain volume : " * string(getDomainVolume(domain)) * " m³")

    timestamps = 0:delta_t:until
    total_frames = length(timestamps)

    # Observable indicating current frame being generated
    # (Used in the record section)
    frame = Observable(1)
    # Observable containing time info
    T = Observable(0.0)

    # Compute all positions for all molecules
    """ pos_hist = [
        mol1: [10.0, 1.4],
        mol2: [13.2, 4.5]
    ]"""
    (pos_hist::Vector{Vector}, speed_hist::Vector{Vector}) = computePositions(mols, domain, delta_t, until)

    # Prepare figure
    fig = Figure(size = (800, 600))
    xlims = (-domain.l_x/2, domain.l_x/2)
    ylims = (-domain.l_y/2, domain.l_y/2)
    zlims = (-domain.l_z/2, domain.l_z/2)

    ax = Axis3(
        fig[1, 1], 
        # perspectiveness = 0.5,
        aspect = (1, 1, 1), 
        title = @lift("t = $($T) s"),
        limits=(xlims, ylims, zlims),
        azimuth = 0.3 * pi
    )

    # Init 1st position for the scatter
    xs, ys, zs = pos_hist[1]
    scatter_plot = scatter!(ax, xs, ys, zs, color=1:20, markersize=10)
    
    record(fig, output_path, 1:framerate:total_frames; framerate = framerate) do f
        frame[] = f
        T[] = timestamps[f]

        # Update scatter with current positions
        xs = []; ys = []; zs = []
        for (i, m) in enumerate(mols) # foreach molecule
            x, y, z = pos_hist[i][f]
            push!(xs, x)
            push!(ys, y)
            push!(zs, z)
        end
        scatter_plot[1] = xs  # xs = [mol1_x, mol2_x, mol3_z] at frame f
        scatter_plot[2] = ys  # y
        scatter_plot[3] = zs  # z
    end

    if (exportToCSV)
        CSV.write("out/output.csv", pos_hist)

        # Export settings
        settings = Dict( 
            "simultation" => Dict(
                "delta_t" => delta_t,
                "until" => until,
                "domain" => [domain.l_x, domain.l_y, domain.l_z],
                "framerate" => framerate
            )
        )

        open("out/output.toml", "w") do io
            TOML.print(io, settings)
        end
    end

    println("Video saved in " * output_path * " ! :)")
end

end # module ISC_MoleculeSim
