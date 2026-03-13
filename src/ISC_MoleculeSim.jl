module ISC_MoleculeSim

using GLMakie, LinearAlgebra, Makie

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

"""Compute all positions for all molecules until provided time."""
function computePositions(mols::Array{Molecule}, domain::Domain, dt::Number, until::Number)::Vector
    pos::Vector{Vector} = []

    # Create new array for molecule if doesn't exist
    for _ in 1:length(mols)
        push!(pos, [])
    end

    for t in 0:dt:until
        # Check Collisions between molecules and walls
        checkWallCollisions(domain, mols)

        # Compute movement
        for (idx, m) in enumerate(mols)
            newPos = computeMovement(m, dt)
            push!(pos[idx], newPos)   # current_m += newPos
        end

        # Check Collisions between molecules
        computeMolsCollisions(mols)
    end

    return pos
end

"""Returns True if m1 overlaps m2."""
function detectMolsCollision(m1::Molecule, m2::Molecule)::Bool
    if (m1 == m2) return false end

    r = m1.radius + m2.radius
    dist = sqrt(sum((m1.position .- m2.position).^2))

    return (dist <= r)
end

"""Generate simulation with provided settings and outputs it to `./out` folder."""
function generateSimulation(domain::Domain, mols::Array{Molecule}, delta_t::Number, until::Number, framerate::Int)::Nothing
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
    pos_hist::Vector{Vector} = computePositions(mols, domain, delta_t, until)

    # Prepare figure
    fig = Figure(size = (800, 600))
    xlims = (-domain.l_x/2, domain.l_x/2)
    ylims = (-domain.l_y/2, domain.l_y/2)
    zlims = (-domain.l_z/2, domain.l_z/2)

    ax = Axis3(
        fig[1, 1], 
        # perspectiveness = 0.5,
        aspect = (1, 1, 1), 
        title = @lift("T = $($T) sec."),
        limits=(xlims, ylims, zlims),
        azimuth = 0.3 * pi
    )

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

    println("Video saved in " * output_path * " ! :)")
end

function main_test()
    """
    Mass of Molecules :
    Helium (He) = 4,002 602 ± 0,000 002 u
    Néon (Ne) = 20,179 7 ± 0,000 6 u
    Azote (N) = 14,006 7 ± 0,000 2 u
    Diazote (N^2) = 2 * N
    Oxygène (O) = 15,999 4 ± 0,000 3 u
    Dioxygène (O^2) = 2 * O  

    Atomic Radius :  
    Helium (He) = 128 pm
    Néon (Ne) = 38 pm
    Diazote (N^2) = 0,315 nm / 2
    Dioxygène (O^2) = 0,292 nm / 2
    (picomètre = 1 * 10^-12 mètres
     nanomètre = 1 * 10^-9 mètres)
    """
    u = 1.660_538_921 * (10^-27) # [kg]
    pico = 10^-12
    nano = 10^-9
    # mHe  = Molecule("He", 4.002_602 * u      , 123 * pico        , [0, -0.0001, -0.0001], [0, 0.003, 0.003])
    # mNe  = Molecule("Ne", 20.179_70 * u      , 38  * pico        , [0, 0.0001, 0.0001], [0, -0.002, -0.002])
    # mN_2 = Molecule("N2", 14.006_70 * 2 * u  , (0.315 / 2) * nano, [0, 0, 0], [0, 0, 0.001])
    # mO_2 = Molecule("O2", 15.999_40 * 2 * u  , (0.292 / 2) * nano, [0, 0, 0], [0, 0, 0.001])

    # mHe  = Molecule("He", 1, 0.5, [0, 0, -4], [0, 0, 2])
    # mNe  = Molecule("Ne", 1, 0.5, [0, 0, 4], [0, 0, -3])

    mHe  = Molecule("He", 1, 0.5, [0, 0, -8], [0, 0, -2])
    mNe  = Molecule("Ne", 1, 0.5, [0, -4, 5], [0, -4, 4])
    mNe2 = Molecule("Ne", 1, 0.5, [0, 0, 4], [0, 0, 4])
    mNe3 = Molecule("Ne", 1, 0.5, [1, 2, 4], [20, 15, 2])

    # Molecules to add to simulation
    molecules::Array{Molecule} = [mHe, mNe, mNe2, mNe3]

    # Domain
    domain::Domain = Domain(20, 20, 20)

    # Time settings
    delta_t::Number = 0.001
    until::Number = 3
    framerate = 30

    generateSimulation(domain, molecules, delta_t, until, framerate)

end

function main_helium()
    # Time settings
    delta_t::Number = 1 *10^-14
    until::Number = 2 *10^-11
    framerate = 30

    # Domain settings
    nano = 10^-9
    size = 10 * nano
    domain::Domain = Domain(size, size, size)

    # Molecules
    num_mols = 400
    molecules::Array{Molecule} = []
    init_speed = 1400 # m/s
    for i in 1:num_mols
        pos::Vector = [
            rand(-1:0.1:1) * domain.l_x/2,
            rand(-1:0.1:1) * domain.l_y/2,
            rand(-1:0.1:1) * domain.l_z/2
        ]

        rand_vect::Vector = randn(3)
        rand_vect = normalize(rand_vect)

        speed::Vector = rand_vect .* init_speed

        mHe = Molecule(
            "He",               # Chemical formula
            6.646 * 10^-27,     # Mass
            1.1 * 10^-10,       # Radius
            pos,                # Position
            speed               # Speed
        )
        push!(molecules, mHe)
    end

    # GENERATE THE AWESOME SIMULATION
    generateSimulation(domain, molecules, delta_t, until, framerate)
end

main_helium()

end # module ISC_MoleculeSim
