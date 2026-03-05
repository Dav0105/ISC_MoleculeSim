module ISC_MoleculeSim

using GLMakie

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

"""Compute next position."""
function computeNextPosition(m::Molecule, dt::Number)::Vector
    m.position = m.position .+ (m.speed * dt)
    return m.position
end

"""Compute all positions for a molecule until provided time."""
function computePositions(m::Molecule, dt::Number, until::Number)::Vector
    pos::Vector = []

    for i in 0:dt:until
        newPos = computeNextPosition(m, dt)
        push!(pos, newPos)
    end

    return pos
end

function detectCollision(m1::Molecule, m2::Molecule)::Bool
    r = m1.radius + m2.radius
    println(m1.radius)
    x1, y1, z1 = m1.position
    x2, y2, z2 = m2.position
    dist = sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)

    return (dist <= r)
end

function main()
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
    mHe  = Molecule("He", 4.002_602 * u      , 123 * pico        , [0, 0, 0], [0.001, 0, 0])
    mNe  = Molecule("Ne", 20.179_70 * u      , 38  * pico        , [0, 0, 0], [0, 0, 0.001])
    mN_2 = Molecule("N2", 14.006_70 * 2 * u  , (0.315 / 2) * nano, [0, 0, 0], [0, 0, 0.001])
    mO_2 = Molecule("O2", 15.999_40 * 2 * u  , (0.292 / 2) * nano, [0, 0, 0], [0, 0, 0.001])

    # Molecules to add to simulation
    molecules::Array = [mHe, mNe]

    # Time settings
    delta_t::Number = 0.001
    until::Number = 3
    framerate = 30

    total_frames = framerate * until
    timestamps = 0:delta_t:until

    # Observable indicating current frame being generated
    # (Used in the record section)
    frame = Observable(1)
    # Observable containing time info
    T = Observable(0.0)

    # Compute all positions for all molecules
    pos_hist::Array{Array} = []
    for m in molecules
        push!(pos_hist, computePositions(m, delta_t, until))
    end

    # Prepare figure
    fig = Figure(size = (800, 600))
    xlims = (-0.0001, 0.0001)
    ylims = (-0.0001, 0.0001)
    zlims = (-0.0001, 0.0001)

    ax = Axis3(
        fig[1, 1], 
        # perspectiveness = 0.5,
        aspect = (1, 1, 1), 
        title = @lift("T = $(round($T, digits=3)) sec."),
        limits=(xlims, ylims, zlims),
        azimuth = 0.3 * pi
    )

    xs, ys, zs = pos_hist[1]
    scatter_plot = scatter!(ax, xs, ys, zs, color=:blue, markersize=10)
    
    record(fig, "out/animation.mp4", 1:total_frames; framerate = framerate) do f
        frame[] = f
        T[] = timestamps[f]

        # Update scatter with current positions
        xs = []; ys = []; zs = []
        for (i, m) in enumerate(molecules)
            x, y, z = pos_hist[i][f]
            push!(xs, x)
            push!(ys, y)
            push!(zs, z)
        end
        scatter_plot[1] = xs  # x
        scatter_plot[2] = ys  # y
        scatter_plot[3] = zs  # z
    end

end

main()

end # module ISC_MoleculeSim
