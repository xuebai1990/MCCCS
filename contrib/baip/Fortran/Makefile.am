# Load the automatically generated dependency file
# include make.depend
rdf.o : global.o
hbond.o : global.o

lib_LTLIBRARIES = libanalysis.la
# ls *.f | awk '{ printf("%s ", $1) }' > tmp and paste the contents of tmp to XXX_SOURCES
libanalysis_la_SOURCES = global.F90 rdf.F90 hbond.F90
libanalysis_la_LDFLAGS = -version-info 0:0:0
