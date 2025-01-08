CFLAGS=-m64 -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -D_BSD_SOURCE -D_POSIX_SOURCE -D_POSIX_C_SOURCE=200809L -D_SVID_SOURCE -D_DARWIN_C_SOURCE -Wall -fno-math-errno -fPIC -std=c99 -lm -O3 -I/usr/include/tirpc 
OFLAGS = -lm -O3 -std=c99 -rdynamic -g -ltirpc
CC = gcc
#CFILES = rockstar.c myutils.c check_syscalls.c fof.c groupies.c subhalo_metric.c potential.c nfw.c jacobi.c fun_times.c interleaving.c universe_time.c hubble.c integrate.c distance.c config_vars.c config.c bounds.c inthash.c io/read_config.c client.c server.c merger.c inet/socket.c inet/rsocket.c inet/address.c io/meta_io.c io/io_internal.c io/io_ascii.c io/stringparse.c io/io_gadget.c io/io_gadgetStar.c io/io_generic.c io/io_art.c io/io_tipsy.c io/io_bgc2.c io/io_util.c io/io_arepo.c io/io_hdf5.c io/getcenter.c io/getrho.c
DIST_FLAGS =
HDF5_FLAGS = -DH5_USE_16_API -lhdf5 -DENABLE_HDF5 -I/opt/local/include -L/opt/local/lib

INCL   = bitarray.h         distance.h      inthash.h    server.h \
	 fof.h           jacobi.h     subhalo_metric.h \
	 bounds.h           fun_times.h     merger.h     universal_constants.h \
	 check_syscalls.h   groupies.h      myutils.h    universe_time.h \
	 client.h           halo.h          nfw.h        version.h \
	 config.h           hubble.h        particle.h \
	 config.template.h  integrate.h     potential.h \
	 config_vars.h      interleaving.h  rockstar.h \
	 io/bgc2.h       io/io_arepo.h  io/io_bgc2.h        io/io_generic.h   io/io_tipsy.h  io/read_config.h \
	 io/getcenter.h  io/io_art.h    io/io_gadget.h      io/io_hdf5.h      io/io_util.h   io/stringparse.h \
	 io/getrho.h     io/io_ascii.h  io/io_gadgetStar.h  io/io_internal.h  io/meta_io.h   io/getv0.h\
	 inet/address.h  inet/rsocket.h inet/socket.h       io/getsigma.h     io/getprof.h io/getMprof.h \
	 util/read_tree.h io/getVprof.h io/getsigmaprof.h   io/io_gadgeth2.h  config_part_mass_ZXY.h\
         ZXY_module/calc_pot_MW.h
# ZXY 2022.02.06: Add io_gadgeth2
# ZXY 2024.03.17: Add particle mass discrimination 
# Note: .o file should also be added below ......
 
EXEC   = modRS

OBJS   = main.o rockstar.o myutils.o check_syscalls.o fof.o groupies.o subhalo_metric.o potential.o nfw.o jacobi.o fun_times.o interleaving.o universe_time.o hubble.o integrate.o distance.o config_vars.o config.o bounds.o inthash.o io/read_config.o client.o server.o merger.o inet/socket.o inet/rsocket.o inet/address.o io/meta_io.o io/io_internal.o io/io_ascii.o io/stringparse.o io/io_gadget.o io/io_gadgetStar.o io/io_generic.o io/io_art.o io/io_tipsy.o io/io_bgc2.o io/io_util.o io/io_arepo.o io/io_hdf5.o io/getcenter.o io/getrho.o io/getv0.o io/getsigma.o io/getprof.o io/getMprof.o io/getVprof.o io/getsigmaprof.o io/io_gadgeth2.o config_part_mass_ZXY.o ZXY_module/calc_pot_MW.o
# ZXY 2022.02.06: Add io_gadgeth2


$(EXEC): $(OBJS)
	 $(CC) $(CFLAGS) $(OBJS) $(LIBS)   -o  $(EXEC) $(OFLAGS)

$(OBJS): $(INCL)
 

#all:
#	@make reg EXTRA_FLAGS="$(OFLAGS)"

#with_hdf5:
#	@make reg EXTRA_FLAGS="$(OFLAGS) $(HDF5_FLAGS)"

#reg:
#	$(CC) $(CFLAGS) main.c $(CFILES) -o modRS1  $(EXTRA_FLAGS)

#lib:
#	$(CC) $(CFLAGS) $(LDFLAGS) $(CFILES) -o librockstar.so $(OFLAGS)

clean:
	rm -f $(OBJS) $(EXEC) *~ io/*~ inet/*~ util/*~ rockstar util/redo_bgc2 util/subhalo_stats

