# Installation directory.
    INSTDIR = $(HIH)

    # List of C++ source files to build.
    # SOP_Main.C registers the operators and handles the DSO-specifics.
    SOURCES = \
	./src/VRAY_InterpolatedGeometry.C


    # Use the highest optimization level.
    OPTIMIZER = -O3

    # Additional include directories.
    INCDIRS = \
        -I/$(HOME)/work/alglib/cpp/include

    # Additional library directories.
    LIBDIRS = \
        -L$(HOME)/work/alglib/lib

    # Additional libraries.
    LIBS = \
        $(HOME)/work/alglib/cpp/lib/interpolation.o \
	$(HOME)/work/alglib/cpp/lib/ap.o \
	$(HOME)/work/alglib/cpp/lib/alglibinternal.o \
	$(HOME)/work/alglib/cpp/lib/alglibmisc.o \
	$(HOME)/work/alglib/cpp/lib/linalg.o \
	$(HOME)/work/alglib/cpp/lib/solvers.o \
	$(HOME)/work/alglib/cpp/lib/optimization.o \
	$(HOME)/work/alglib/cpp/lib/specialfunctions.o \
	$(HOME)/work/alglib/cpp/lib/integration.o

    # Set the plugin library name.
    DSONAME = VRAY_InterpolatedGeometry.so

    # Include the GNU Makefile.
    include $(HFS)/toolkit/makefiles/Makefile.gnu

