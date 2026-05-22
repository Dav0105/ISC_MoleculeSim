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
    
    mass::Float64        # kg
    radius::Float64      # m

    position::Vector{Float64}
    speed::Vector{Float64}       # m/s

    pos_hist::Array{Vector{Float64}}
    speed_hist::Array{Vector{Float64}}
end

"""Represents the volume where the molecules will be."""
struct Domain
    l_x::Float64 # m
    l_y::Float64 # m
    l_z::Float64 # m
end

"""Returns the volume of the provided domain in m³."""
function getDomainVolume(d::Domain)::Number
    return d.l_x * d.l_y * d.l_z
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
function computePositions(mols::Array{Molecule}, domain::Domain, dt::Number, until::Number, g::Number)

    for t in 0:dt:until
        # Check Collisions between molecules and walls
        checkWallCollisions(domain, mols)

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

#region Graphs generation

"""Displays a Makie graph containing the distributions of positions 
in Z of the molecules, at the end of the simulation."""
function generateGravityProbaGraph(mols::Array{Molecule}, filename="./out/z_position_dist.png")
    last_pos_hist = []
    for mol in mols
        append!(last_pos_hist, last(mol.pos_hist)[3])
    end

    f = Figure()
    ax = Axis(f[1, 1],
        title = "Distribution of molecules positions in Z at the end of the simulation",
        xlabel = "Position in Z",
        ylabel = "Number of molecules",
    )
    hist!(ax, last_pos_hist)
    save(filename, f)           # Save figure
    display(f)
    return
end

#endregion

"""Generate simulation with provided settings and outputs it to `./out` folder."""
function generateSimulation(domain::Domain, mols::Array{Molecule}, delta_t::Number, until::Number, framerate::Int; framestep::Int= 30, exportToCSV::Bool = true, output_path::String = "./out/animation", g::Number=-9.81)    
    println("Domain volume : " * string(getDomainVolume(domain)) * " m³")

    timestamps = 0:delta_t:until
    total_frames = length(timestamps)

    # Observable indicating current frame being generated
    # (Used in the record section)
    frame = Observable(1)
    # Observable containing time info
    T = Observable(0.0)

    @time "Time to generate positions" begin
        # Compute all positions for all molecules
        computePositions(mols, domain, delta_t, until, g)
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
                "domain" => [domain.l_x, domain.l_y, domain.l_z],
                "framerate" => framerate,
                "g" => g
            )
        )
        open(output_path * ".toml", "w") do io
            TOML.print(io, settings)
        end
    end


    println("Video saved in " * output_path * " ! :)")

    generateGravityProbaGraph(mols, output_path * "_z_positions_dist.png")

    display(fig)

end

end # module ISC_MoleculeSim
