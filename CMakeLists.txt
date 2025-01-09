cmake_minimum_required(VERSION 3.30)
project(cosmoAnalysis C)

set(CMAKE_C_STANDARD 11)

include_directories(commands)
include_directories(examples)
include_directories(inet)
include_directories(io)
include_directories(scripts)
include_directories(util)
include_directories(ZXY_module)

## Automatically collect all .h and .c files
#file(GLOB PROJECT_HEADERS "RS/*.h")
#file(GLOB PROJECT_SOURCES "RS/*.c")
#
## Add them to your executable or library
#add_executable(cosmoAnalysis ${PROJECT_HEADERS} ${PROJECT_SOURCES}
#        ../Gadget-2/cosmoRS/scripts/runsingle.c)


#add_executable(cosmoAnalysis
#        main.c)

# Updated source and header files
add_executable(cosmoAnalysis
    commands/count_CDM760.c
    commands/count_presim.c
    commands/mvfile.c
    commands/getbdmass_CDM_SIDM.c
    commands/getDMpos_CDM.c
    commands/count_SIDM100.c
    commands/getDMprof_CDM.c
    commands/count_SIDM.c
    commands/getDMpos_SIDM.c
    commands/getDMprof_SIDM.c
    examples/load_bgc2.c
    examples/bgc2_to_ascii_particles.c
    examples/calc_potentials.c
    examples/load_full_particles.c
    examples/calc_bgc2_shapes.c
    examples/load_full_particles.h
    examples/load_bgc2.h
    inet/rsocket.c
    inet/socket.c
    inet/address.c
    inet/address.h
    inet/rsocket.h
    inet/socket.h
    io/extra_bgc2.c
    io/io_generic.c
    io/io_arepo.c
    io/meta_io.c
    io/io_gadgetStar.c
    io/io_gadget.c
    io/io_bgc2.c
    io/getv0.c
    io/getsigmaprof.c
    io/io_gadgeth2.c
    io/getrho.c
    io/getsigma.c
    io/getvr.c
    io/io_hdf5.c
    io/stringparse.c
    io/strtonum.c
    io/io_art.c
    io/getVprof.c
    io/getMprof.c
    io/io_ascii.c
    io/getsigmaprofZXY.c
    io/getcenter.c
    io/io_util.c
    io/getprof.c
    io/io_internal.c
    io/io_tipsy.c
    io/read_config.c
    io/stringparse.h
    io/io_ascii.h
    io/bgc2.h
    io/getMprof.h
    io/getVprof.h
    io/io_art.h
    io/io_util.h
    io/getprof.h
    io/getcenter.h
    io/io_tipsy.h
    io/read_config.h
    io/io_internal.h
    io/io_arepo.h
    io/io_generic.h
    io/io_gadgetStar.h
    io/io_gadget.h
    io/meta_io.h
    io/getsigma.h
    io/getrho.h
    io/io_gadgeth2.h
    io/getsigmaprof.h
    io/getv0.h
    io/io_bgc2.h
    io/io_hdf5.h
    scripts/runsingle.c
    util/redo_bgc2.c
    util/bgc2_to_ascii.c
    util/find_parents.c
    util/subhalo_stats.c
    util/read_tree.h
    ZXY_module/gravtree.c
    ZXY_module/calc_pot_MW.c
    ZXY_module/calc_pot_MW.h
    potential.c
    merger.c
    fof.c
    nfw.c
    bounds.c
    load_balance.c
    jacobi.c
    subhalo_metric.c
    fun_times.c
    config_vars.c
    singleprop.c
    interleaving.c
    myutils.c
    properties.c
    rockstar.c
    server.c
    client.c
    inthash.c
    main.c
    integrate.c
    config.c
    parents.c
    distance.c
    hubble.c
    config_part_mass_ZXY.c
    starprop.c
    fast3tree.c
    universe_time.c
    groupies.c
    check_syscalls.c
    bitarray.h
    inthash.h
    client.h
    version.h
    hubble.h
    distance.h
    config.h
    integrate.h
    halo.h
    config_part_mass_ZXY.h
    check_syscalls.h
    groupies.h
    universe_time.h
    nfw.h
    fof.h
    merger.h
    potential.h
    jacobi.h
    particle.h
    bounds.h
    myutils.h
    config_vars.h
    fun_times.h
    interleaving.h
    subhalo_metric.h
    server.h
    universal_constants.h
    bkp_singleprop.h
    config.template.h
    rockstar.h)