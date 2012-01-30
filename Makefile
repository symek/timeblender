# Installation directory.
    INSTDIR = $(HIH)

    # List of C++ source files to build.
    # SOP_Main.C registers the operators and handles the DSO-specifics.
    SOURCES = ./src/VRAY_TimeBlender.C ./src/TB_PointMatch.C ./src/TB_GeoInterpolants.C


    # Use the highest optimization level.
    OPTIMIZER = -O3

    # Additional include directories.
    #INCDIRS = \
    #   -I/$(HOME)/

    # Additional library directories.
    # LIBDIRS = \
    #    -L$(HOME)/
   
    # Additional libraries.
    #LIBS = ./src/TB_PointMatch.o
	
    # Set the plugin library name.
    DSONAME = VRAY_TimeBlender.so

    # Include the GNU Makefile.
    include $(HFS)/toolkit/makefiles/Makefile.gnu
