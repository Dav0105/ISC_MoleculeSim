using LinearAlgebra
# using ISC_MoleculeSim

include("ISC_MoleculeSim.jl")

function generate_random_positions(domain::ISC_MoleculeSim.Domain)
    return [
        rand(-1:0.001:1) * (domain.lims_x[2] - domain.lims_x[1]) / 2 + (domain.lims_x[1] + domain.lims_x[2]) / 2,
        rand(-1:0.001:1) * (domain.lims_y[2] - domain.lims_y[1]) / 2 + (domain.lims_y[1] + domain.lims_y[2]) / 2,
        rand(-1:0.001:1) * (domain.lims_z[2] - domain.lims_z[1]) / 2 + (domain.lims_z[1] + domain.lims_z[2]) / 2
    ]
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

    mHe  = ISC_MoleculeSim.Molecule("He", 1, 0.5, [0, 0, -8], [0, 0, -2], [], [])
    mNe  = ISC_MoleculeSim.Molecule("Ne", 1, 0.5, [0, -4, 5], [0, -4, 4], [], [])
    mNe2 = ISC_MoleculeSim.Molecule("Ne", 1, 0.5, [0, 0, 4], [0, 0, 4], [], [])
    mNe3 = ISC_MoleculeSim.Molecule("Ne", 1, 0.5, [1, 2, 4], [20, 15, 2], [], [])

    # Molecules to add to simulation
    molecules::Array{ISC_MoleculeSim.Molecule} = [mHe, mNe, mNe2, mNe3]

    # Domain
    domain::ISC_MoleculeSim.Domain = ISC_MoleculeSim.Domain(
        (-10, 10),
        (-10, 10),
        (-10, 10)
    )

    # Time settings
    delta_t::Number = 0.0001
    until::Number = 3
    framerate = 30

    ISC_MoleculeSim.generateSimulation(domain, molecules, delta_t, until, framerate)

end

function main_helium()
    # Time settings
    delta_t::Number = 1 *10^-14
    until::Number = 2 *10^-11
    framerate = 30

    # Domain settings
    nano = 10^-9
    size = 10 * nano
    domain::ISC_MoleculeSim.Domain = ISC_MoleculeSim.Domain(
        (-size / 2, size / 2), 
        (-size / 2, size / 2), 
        (-size / 2, size / 2)
    )

    # Molecules
    num_mols = 400
    molecules::Array{ISC_MoleculeSim.Molecule} = []
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

        mHe = ISC_MoleculeSim.Molecule(
            "He",               # Chemical formula
            6.646 * 10^-27,     # Mass
            1.1 * 10^-10,       # Radius
            pos,                # Position
            speed,              # Speed
            [],                 # pos_hist
            []                  # speed_hist
        )
        push!(molecules, mHe)
    end

    # GENERATE THE AWESOME SIMULATION
    ISC_MoleculeSim.generateSimulation(domain, molecules, delta_t, until, framerate)
end

function main_helium_2()
    # Time settings
    delta_t::Number = 1 *10^-14
    until::Number = 5 * 10^-11
    framerate = 30

    # Domain settings
    size = 2 * 10^-8
    domain::ISC_MoleculeSim.Domain = ISC_MoleculeSim.Domain((-size / 2, size / 2), (-size / 2, size / 2), (-size / 2, size / 2))

    # g
    g = -9.81 * 10^14

    # Molecules
    molecules::Array{ISC_MoleculeSim.Molecule} = []

    ## Helium
    num_mols_helium = 400
    init_speed_helium = 789.45 # m/s
    for i in 1:num_mols_helium
        pos::Vector = generate_random_positions(domain)

        # Generate random speeds
        rand_vect::Vector = randn(3)
        rand_vect = normalize(rand_vect)
        speed::Vector = rand_vect .* init_speed_helium

        mHe = ISC_MoleculeSim.Molecule(
            "He",               # Chemical formula
            6.646 * 10^-27,     # Mass
            1.1 * 10^-10,       # Radius
            pos,                # Position
            speed,              # Speed
            [],                 # pos_hist
            []                  # speed_hist
        )
        push!(molecules, mHe)
    end

    ## Argon
    num_mols_argon = 200
    init_speed_argon = 249.88 # m/s
    for i in 1:num_mols_argon
        pos::Vector = generate_random_positions(domain)

        # Generate random speeds
        rand_vect::Vector = randn(3)
        rand_vect = normalize(rand_vect)
        speed::Vector = rand_vect .* init_speed_argon

        mAr = ISC_MoleculeSim.Molecule(
            "Ar",               # Chemical formula
            6.634 *10^-26,      # Mass
            1.88 * 10^-10,      # Radius
            pos,                # Position
            speed,              # Speed
            [],                 # pos_hist
            []                  # speed_hist
        )
        push!(molecules, mAr)
    end

    # GENERATE THE AWESOME SIMULATION
    ISC_MoleculeSim.generateSimulation(domain, molecules, delta_t, until, framerate, framestep=30, g=g)
end

function main_helium_2_lite()
    # Time settings
    delta_t::Number = 1 *10^-14
    until::Number = 10 * 10^-11
    framerate = 30

    # Domain settings
    size = 0.25 * 10^-8
    domain::ISC_MoleculeSim.Domain = ISC_MoleculeSim.Domain((-size / 2, size / 2), (-size / 2, size / 2), (-size / 2, size / 2))

    # g
    g = -9.81 * 10^13

    # Molecules
    molecules::Array{ISC_MoleculeSim.Molecule} = []

    ## Helium
    num_mols_helium = 200
    init_speed_helium = 789.45 # m/s
    for i in 1:num_mols_helium
        pos::Vector = generate_random_positions(domain)

        # Generate random speeds
        rand_vect::Vector = randn(3)
        rand_vect = normalize(rand_vect)
        speed::Vector = rand_vect .* init_speed_helium

        mHe = ISC_MoleculeSim.Molecule(
            "He",               # Chemical formula
            6.646 * 10^-27,     # Mass
            1.1 * 10^-10,       # Radius
            pos,                # Position
            speed,              # Speed
            [],                 # pos_hist
            []                  # speed_hist
        )
        push!(molecules, mHe)
    end

    ## Argon
    num_mols_argon = 100
    init_speed_argon = 249.88 # m/s
    for i in 1:num_mols_argon
        pos::Vector = generate_random_positions(domain)

        # Generate random speeds
        rand_vect::Vector = randn(3)
        rand_vect = normalize(rand_vect)
        speed::Vector = rand_vect .* init_speed_argon

        mAr = ISC_MoleculeSim.Molecule(
            "Ar",               # Chemical formula
            6.634 *10^-26,      # Mass
            1.88 * 10^-10,      # Radius
            pos,                # Position
            speed,              # Speed
            [],                 # pos_hist
            []                  # speed_hist
        )
        push!(molecules, mAr)
    end

    # GENERATE THE AWESOME SIMULATION
    ISC_MoleculeSim.generateSimulation(domain, molecules, delta_t, until, framerate, framestep=30, g=g)
end

function main_helium_small_init_pos()
    # Time settings
    delta_t::Number = 1 *10^-14
    until::Number = 10 * 10^-11
    framerate = 30

    # Domain settings
    domain::ISC_MoleculeSim.Domain = ISC_MoleculeSim.Domain(
        (-5 * 10^-9, 5 * 10^-9), 
        (-5 * 10^-9, 5 * 10^-9), 
        (-5 * 10^-9, 5 * 10^-9)
    )
    spawn_domain::ISC_MoleculeSim.Domain = ISC_MoleculeSim.Domain(
        (-5 * 10^-9 / 8, 5 * 10^-9 / 8), 
        (-5 * 10^-9 / 8, 5 * 10^-9 / 8), 
        (-5 * 10^-9 / 8, 5 * 10^-9 / 8)
    )

    # g
    # g = -9.81 * 10^13
    g = 0.0

    # Molecules
    molecules::Array{ISC_MoleculeSim.Molecule} = []

    ## Helium
    num_mols_helium = 400
    init_speed_helium = 1400.0 # m/s
    for i in 1:num_mols_helium
        pos::Vector = generate_random_positions(spawn_domain)

        # Generate random speeds
        rand_vect::Vector = randn(3)
        rand_vect = normalize(rand_vect)
        speed::Vector = rand_vect .* init_speed_helium

        mHe = ISC_MoleculeSim.Molecule(
            "He",               # Chemical formula
            6.646 * 10^-27,     # Mass
            1.1 * 10^-10,       # Radius
            pos,                # Position
            speed,              # Speed
            [],                 # pos_hist
            []                  # speed_hist
        )
        push!(molecules, mHe)
    end

    # GENERATE THE AWESOME SIMULATION
    ISC_MoleculeSim.generateSimulation(
        domain, molecules, delta_t, until, framerate, 
        framestep=30, g=g, probability_bins=(20, 10, 10, 200)
    )
end

function main_helium_small_expanding_domain()
    # Time settings
    delta_t::Number = 1 *10^-14
    until::Number = 10 * 10^-11
    framerate = 30

    # Domain settings
    domain::ISC_MoleculeSim.Domain = ISC_MoleculeSim.Domain(
        (-5 * 10^-9, 5 * 10^-9), 
        (-5 * 10^-9, 5 * 10^-9), 
        (-5 * 10^-9, 5 * 10^-9)
    )
    second_domain::ISC_MoleculeSim.Domain = ISC_MoleculeSim.Domain(
        (-5 * 10^-9, 15 * 10^-9), 
        (-5 * 10^-9, 5 * 10^-9), 
        (-5 * 10^-9, 5 * 10^-9)
    )
    spawn_domain::ISC_MoleculeSim.Domain = ISC_MoleculeSim.Domain(
        (-5 * 10^-9 / 8, 5 * 10^-9 / 8), 
        (-5 * 10^-9 / 8, 5 * 10^-9 / 8), 
        (-5 * 10^-9 / 8, 5 * 10^-9 / 8)
    )

    # g
    # g = -9.81 * 10^13
    g = 0.0

    # Molecules
    molecules::Array{ISC_MoleculeSim.Molecule} = []

    ## Helium
    num_mols_helium = 400
    init_speed_helium = 1400.0 # m/s
    for i in 1:num_mols_helium
        pos::Vector = generate_random_positions(spawn_domain)

        # Generate random speeds
        rand_vect::Vector = randn(3)
        rand_vect = normalize(rand_vect)
        speed::Vector = rand_vect .* init_speed_helium

        mHe = ISC_MoleculeSim.Molecule(
            "He",               # Chemical formula
            6.646 * 10^-27,     # Mass
            1.1 * 10^-10,       # Radius
            pos,                # Position
            speed,              # Speed
            [],                 # pos_hist
            []                  # speed_hist
        )
        push!(molecules, mHe)
    end

    # GENERATE THE AWESOME SIMULATION
    ISC_MoleculeSim.generateSimulation(
        domain, molecules, delta_t, until, framerate, 
        framestep=30, g=g, second_domain=second_domain, probability_bins=(20, 10, 10, 200)
    )
end

function main_helium_temp()
    # Time settings
    delta_t::Number = 1 *10^-14
    until::Number = 10 * 10^-11
    framerate = 30

    # Domain settings
    domain::ISC_MoleculeSim.Domain = ISC_MoleculeSim.Domain(
        (-2 * 10^-9, 2 * 10^-9), 
        (-2 * 10^-9, 2 * 10^-9), 
        (-5 * 10^-9, 5 * 10^-9)
    )

    # Temperatures (z+, z-)
    temperatures = (700.0, 300.0)

    # g
    # g = -9.81 * 10^13
    g = 0.0

    # Molecules
    molecules::Array{ISC_MoleculeSim.Molecule} = []

    ## Helium
    num_mols_helium = 500
    init_speed_helium = 1400.0 # m/s
    for i in 1:num_mols_helium
        pos::Vector = generate_random_positions(domain)

        # Generate random speeds
        rand_vect::Vector = randn(3)
        rand_vect = normalize(rand_vect)
        speed::Vector = rand_vect .* init_speed_helium

        mHe = ISC_MoleculeSim.Molecule(
            "He",               # Chemical formula
            6.646 * 10^-27,     # Mass
            1.1 * 10^-10,       # Radius
            pos,                # Position
            speed,              # Speed
            [],                 # pos_hist
            []                  # speed_hist
        )
        push!(molecules, mHe)
    end

    # GENERATE THE AWESOME SIMULATION
    ISC_MoleculeSim.generateSimulation(
        domain, molecules, delta_t, until, framerate, 
        framestep=30, g=g, temperatures=temperatures, probability_bins=(20, 10, 10, 200)
    )
end

# main_test()
# main_helium()
# main_helium_2()
# main_helium_2_lite()
main_helium_small_init_pos()
# main_helium_small_expanding_domain()
# main_helium_temp()