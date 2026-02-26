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
function computeNextPosition(m::Molecule, dt::Number)
    m.position = m.position .+ (m.speed * dt)
end

function computePositions(m::Molecule, dt::Number, until::Number)::Vector
    pos::Vector = []

    for i in 0:dt:until
        newPos = computeNextPosition(m, dt)
        push!(pos, newPos)
    end

    return pos
end

function main()
    """
    Mass of Molecules :
    Helium (He) = 4,002 602 ± 0,000 002 u
    Néon (Ne) = 20,179 7 ± 0,000 6 u
    Azote (N) = 14,006 7 ± 0,000 2 uMolecule
    Diazote (N^2) = 2 * N
    Oxygène (O) = 15,999 4 ± 0,000 3 u
    Dioxygène (O^2) = 2 * O
    """
    u = 1.660_538_921 * (10^-27) # [kg]
    mHe  = Molecule("He", 4.002_602 * u      , 0, [0, 0, 0], [0, 0.001, 0.000])
    mNe  = Molecule("Ne", 20.179_70 * u      , 0, [0, 0, 0], [0, 0, 0.001])
    mN_2 = Molecule("N2", 14.006_70 * 2 * u  , 0, [0, 0, 0], [0, 0, 0.001])
    mO_2 = Molecule("O2", 15.999_40 * 2 * u  , 0, [0, 0, 0], [0, 0, 0.001])

    # Time settings
    delta_t::Number = 0.001
    until::Number = 3
    framerate = 30

    # Observable indicating current frame being generated
    # (Used in the record section)
    frame = Observable(1)

    # Compute all positions
    pos_hist = computePositions(mHe, delta_t, 10)

    # Prepare 
    fig = Figure(size = (800, 600))
    xlims = (-0.001, 0.001)
    ylims = (-0.001, 0.001)
    zlims = (-0.001, 0.001)

    Axis3(
        fig[1, 1], 
        # perspectiveness = 0.5,
        aspect = (1, 1, 1), 
        title = @lift(),
        limits=(xlims, ylims, zlims),
        azimuth = 0.3 * pi
    )

    @lift(scatter!(pos_hist[$frame][1], pos_hist[$frame][2], pos_hist[$frame][3]))
    
    total_frames = framerate * until
    timestamps = 0:delta_t:until
    
    record(fig, "animation.mp4", 1:total_frames; framerate = framerate) do f
        frame[] = f
    end

end

main()

end # module ISC_MoleculeSim
