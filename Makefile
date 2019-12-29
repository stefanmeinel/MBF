CXX           = g++
CXXFLAGS      = -O3 -Wall
INCPATH       =
LIBS          = -lgsl -lgslcblas -lm
LINK          = g++
LFLAGS        =
DEL_FILE      = rm -f

TARGET        = MBF

SOURCES       = main.cpp \
		settingsmap.cpp \
		binary_io.cpp \
		parse_model.cpp \
		multi_alt_exp_Asqr_BC_model.cpp \
		multi_alt_exp_Asqr_expE_BC_model.cpp \
		multi_alt_exp_Asqr_expE_model.cpp \
		multi_alt_exp_Asqr_expE_vec_BC_model.cpp \
		multi_alt_exp_Asqr_expE_vec_model.cpp \
		multi_alt_exp_Asqr_model.cpp \
		multi_alt_exp_Asqr_vec_BC_model.cpp \
		multi_alt_exp_Asqr_vec_model.cpp \
		multi_alt_exp_BC_model.cpp \
		multi_alt_exp_expE_BC_model.cpp \
		multi_alt_exp_expE_mat_model.cpp \
		multi_alt_exp_expE_model.cpp \
		multi_alt_exp_expE_nonsym_mat_model.cpp \
		multi_alt_exp_expE_vec_BC_model.cpp \
		multi_alt_exp_expE_vec_model.cpp \
		multi_alt_exp_mat_model.cpp \
		multi_alt_exp_model.cpp \
		multi_alt_exp_nonsym_mat_model.cpp \
		multi_alt_exp_vec_BC_model.cpp \
		multi_alt_exp_vec_model.cpp \
		multi_exp_Asqr_BC_model.cpp \
		multi_exp_Asqr_expE_BC_model.cpp \
		multi_exp_Asqr_expE_model.cpp \
		multi_exp_Asqr_expE_vec_BC_model.cpp \
		multi_exp_Asqr_expE_vec_model.cpp \
		multi_exp_Asqr_model.cpp \
		multi_exp_Asqr_vec_BC_model.cpp \
		multi_exp_Asqr_vec_model.cpp \
		multi_exp_BC_model.cpp \
		multi_exp_expE_BC_model.cpp \
		multi_exp_expE_mat_model.cpp \
		multi_exp_expE_mat_II_model.cpp \
		multi_exp_expE_mat_upper_model.cpp \
		multi_exp_expE_mat_II_upper_model.cpp \
		multi_exp_expE_model.cpp \
		multi_exp_expE_nonsym_mat_model.cpp \
		multi_exp_expE_vec_BC_model.cpp \
		multi_exp_expE_vec_model.cpp \
		multi_exp_mat_model.cpp \
		multi_exp_mat_II_model.cpp \
		multi_exp_mat_upper_model.cpp \
		multi_exp_mat_II_upper_model.cpp \
		multi_exp_model.cpp \
		multi_exp_nonsym_mat_model.cpp \
		multi_exp_vec_BC_model.cpp \
		multi_exp_vec_model.cpp \
		multi_exp_const_model.cpp \
		multi_exp_expE_const_model.cpp \
		multi_exp_Asqr_const_model.cpp \
		multi_exp_Asqr_expE_const_model.cpp \
		multi_alt_exp_const_model.cpp \
		multi_alt_exp_expE_const_model.cpp \
		multi_alt_exp_Asqr_const_model.cpp \
		multi_alt_exp_Asqr_expE_const_model.cpp \
		multi_exp_vec_const_model.cpp \
		multi_exp_expE_vec_const_model.cpp \
		multi_exp_Asqr_vec_const_model.cpp \
		multi_exp_Asqr_expE_vec_const_model.cpp \
		multi_alt_exp_vec_const_model.cpp \
		multi_alt_exp_expE_vec_const_model.cpp \
		multi_alt_exp_Asqr_vec_const_model.cpp \
		multi_alt_exp_Asqr_expE_vec_const_model.cpp \
		multi_exp_BC_const_model.cpp \
		multi_exp_expE_BC_const_model.cpp \
		multi_exp_Asqr_BC_const_model.cpp \
		multi_exp_Asqr_expE_BC_const_model.cpp \
		multi_alt_exp_BC_const_model.cpp \
		multi_alt_exp_expE_BC_const_model.cpp \
		multi_alt_exp_Asqr_BC_const_model.cpp \
		multi_alt_exp_Asqr_expE_BC_const_model.cpp \
		multi_exp_vec_BC_const_model.cpp \
		multi_exp_expE_vec_BC_const_model.cpp \
		multi_exp_Asqr_vec_BC_const_model.cpp \
		multi_exp_Asqr_expE_vec_BC_const_model.cpp \
		multi_alt_exp_vec_BC_const_model.cpp \
		multi_alt_exp_expE_vec_BC_const_model.cpp \
		multi_alt_exp_Asqr_vec_BC_const_model.cpp \
		multi_alt_exp_Asqr_expE_vec_BC_const_model.cpp \
		threept_multi_alt_exp_expE_model.cpp \
		threept_multi_alt_exp_expE_vec_model.cpp \
		threept_multi_alt_exp_model.cpp \
		threept_multi_alt_exp_vec_model.cpp \
		threept_multi_exp_expE_model.cpp \
		threept_multi_exp_expE_vec_model.cpp \
		threept_multi_exp_model.cpp \
		threept_multi_exp_vec_model.cpp \
		gaussian_prior.cpp \
		fitter.cpp \
		parser.cpp \
		chisqr_extra_term.cpp


OBJECTS       = main.o \
		settingsmap.o \
		binary_io.o \
		parse_model.o \
		multi_alt_exp_Asqr_BC_model.o \
		multi_alt_exp_Asqr_expE_BC_model.o \
		multi_alt_exp_Asqr_expE_model.o \
		multi_alt_exp_Asqr_expE_vec_BC_model.o \
		multi_alt_exp_Asqr_expE_vec_model.o \
		multi_alt_exp_Asqr_model.o \
		multi_alt_exp_Asqr_vec_BC_model.o \
		multi_alt_exp_Asqr_vec_model.o \
		multi_alt_exp_BC_model.o \
		multi_alt_exp_expE_BC_model.o \
		multi_alt_exp_expE_mat_model.o \
		multi_alt_exp_expE_model.o \
		multi_alt_exp_expE_nonsym_mat_model.o \
		multi_alt_exp_expE_vec_BC_model.o \
		multi_alt_exp_expE_vec_model.o \
		multi_alt_exp_mat_model.o \
		multi_alt_exp_model.o \
		multi_alt_exp_nonsym_mat_model.o \
		multi_alt_exp_vec_BC_model.o \
		multi_alt_exp_vec_model.o \
		multi_exp_Asqr_BC_model.o \
		multi_exp_Asqr_expE_BC_model.o \
		multi_exp_Asqr_expE_model.o \
		multi_exp_Asqr_expE_vec_BC_model.o \
		multi_exp_Asqr_expE_vec_model.o \
		multi_exp_Asqr_model.o \
		multi_exp_Asqr_vec_BC_model.o \
		multi_exp_Asqr_vec_model.o \
		multi_exp_BC_model.o \
		multi_exp_expE_BC_model.o \
		multi_exp_expE_mat_model.o \
		multi_exp_expE_mat_II_model.o \
		multi_exp_expE_mat_upper_model.o \
		multi_exp_expE_mat_II_upper_model.o \
		multi_exp_expE_model.o \
		multi_exp_expE_nonsym_mat_model.o \
		multi_exp_expE_vec_BC_model.o \
		multi_exp_expE_vec_model.o \
		multi_exp_mat_model.o \
		multi_exp_mat_II_model.o \
		multi_exp_mat_upper_model.o \
		multi_exp_mat_II_upper_model.o \
		multi_exp_model.o \
		multi_exp_nonsym_mat_model.o \
		multi_exp_vec_BC_model.o \
		multi_exp_vec_model.o \
		multi_exp_const_model.o \
		multi_exp_expE_const_model.o \
		multi_exp_Asqr_const_model.o \
		multi_exp_Asqr_expE_const_model.o \
		multi_alt_exp_const_model.o \
		multi_alt_exp_expE_const_model.o \
		multi_alt_exp_Asqr_const_model.o \
		multi_alt_exp_Asqr_expE_const_model.o \
		multi_exp_vec_const_model.o \
		multi_exp_expE_vec_const_model.o \
		multi_exp_Asqr_vec_const_model.o \
		multi_exp_Asqr_expE_vec_const_model.o \
		multi_alt_exp_vec_const_model.o \
		multi_alt_exp_expE_vec_const_model.o \
		multi_alt_exp_Asqr_vec_const_model.o \
		multi_alt_exp_Asqr_expE_vec_const_model.o \
		multi_exp_BC_const_model.o \
		multi_exp_expE_BC_const_model.o \
		multi_exp_Asqr_BC_const_model.o \
		multi_exp_Asqr_expE_BC_const_model.o \
		multi_alt_exp_BC_const_model.o \
		multi_alt_exp_expE_BC_const_model.o \
		multi_alt_exp_Asqr_BC_const_model.o \
		multi_alt_exp_Asqr_expE_BC_const_model.o \
		multi_exp_vec_BC_const_model.o \
		multi_exp_expE_vec_BC_const_model.o \
		multi_exp_Asqr_vec_BC_const_model.o \
		multi_exp_Asqr_expE_vec_BC_const_model.o \
		multi_alt_exp_vec_BC_const_model.o \
		multi_alt_exp_expE_vec_BC_const_model.o \
		multi_alt_exp_Asqr_vec_BC_const_model.o \
		multi_alt_exp_Asqr_expE_vec_BC_const_model.o \
		threept_multi_alt_exp_expE_model.o \
		threept_multi_alt_exp_expE_vec_model.o \
		threept_multi_alt_exp_model.o \
		threept_multi_alt_exp_vec_model.o \
		threept_multi_exp_expE_model.o \
		threept_multi_exp_expE_vec_model.o \
		threept_multi_exp_model.o \
		threept_multi_exp_vec_model.o \
		gaussian_prior.o \
		fitter.o \
		parser.o \
		chisqr_extra_term.o


$(TARGET):  $(OBJECTS)
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(LIBS)


main.o: 	main.cpp \
		fitter.h \
		abstract_model.h \
		abstract_prior.h \
		gaussian_prior.h \
		parse_model.h \
		parser.h \
		chisqr_extra_term.h \
		binary_io.h \
		multi_alt_exp_Asqr_BC_model.h \
		multi_alt_exp_Asqr_expE_BC_model.h \
		multi_alt_exp_Asqr_expE_model.h \
		multi_alt_exp_Asqr_expE_vec_BC_model.h \
		multi_alt_exp_Asqr_expE_vec_model.h \
		multi_alt_exp_Asqr_model.h \
		multi_alt_exp_Asqr_vec_BC_model.h \
		multi_alt_exp_Asqr_vec_model.h \
		multi_alt_exp_BC_model.h \
		multi_alt_exp_expE_BC_model.h \
		multi_alt_exp_expE_mat_model.h \
		multi_alt_exp_expE_model.h \
		multi_alt_exp_expE_nonsym_mat_model.h \
		multi_alt_exp_expE_vec_BC_model.h \
		multi_alt_exp_expE_vec_model.h \
		multi_alt_exp_mat_model.h \
		multi_alt_exp_model.h \
		multi_alt_exp_nonsym_mat_model.h \
		multi_alt_exp_vec_BC_model.h \
		multi_alt_exp_vec_model.h \
		multi_exp_Asqr_BC_model.h \
		multi_exp_Asqr_expE_BC_model.h \
		multi_exp_Asqr_expE_model.h \
		multi_exp_Asqr_expE_vec_BC_model.h \
		multi_exp_Asqr_expE_vec_model.h \
		multi_exp_Asqr_model.h \
		multi_exp_Asqr_vec_BC_model.h \
		multi_exp_Asqr_vec_model.h \
		multi_exp_BC_model.h \
		multi_exp_expE_BC_model.h \
		multi_exp_expE_mat_model.h \
		multi_exp_expE_mat_II_model.h \
		multi_exp_expE_mat_upper_model.h \
		multi_exp_expE_mat_II_upper_model.h \
		multi_exp_expE_model.h \
		multi_exp_expE_nonsym_mat_model.h \
		multi_exp_expE_vec_BC_model.h \
		multi_exp_expE_vec_model.h \
		multi_exp_mat_model.h \
		multi_exp_mat_II_model.h \
		multi_exp_mat_upper_model.h \
		multi_exp_mat_II_upper_model.h \
		multi_exp_model.h \
		multi_exp_nonsym_mat_model.h \
		multi_exp_vec_BC_model.h \
		multi_exp_vec_model.h \
		multi_exp_const_model.h \
		multi_exp_expE_const_model.h \
		multi_exp_Asqr_const_model.h \
		multi_exp_Asqr_expE_const_model.h \
		multi_alt_exp_const_model.h \
		multi_alt_exp_expE_const_model.h \
		multi_alt_exp_Asqr_const_model.h \
		multi_alt_exp_Asqr_expE_const_model.h \
		multi_exp_vec_const_model.h \
		multi_exp_expE_vec_const_model.h \
		multi_exp_Asqr_vec_const_model.h \
		multi_exp_Asqr_expE_vec_const_model.h \
		multi_alt_exp_vec_const_model.h \
		multi_alt_exp_expE_vec_const_model.h \
		multi_alt_exp_Asqr_vec_const_model.h \
		multi_alt_exp_Asqr_expE_vec_const_model.h \
		multi_exp_BC_const_model.h \
		multi_exp_expE_BC_const_model.h \
		multi_exp_Asqr_BC_const_model.h \
		multi_exp_Asqr_expE_BC_const_model.h \
		multi_alt_exp_BC_const_model.h \
		multi_alt_exp_expE_BC_const_model.h \
		multi_alt_exp_Asqr_BC_const_model.h \
		multi_alt_exp_Asqr_expE_BC_const_model.h \
		multi_exp_vec_BC_const_model.h \
		multi_exp_expE_vec_BC_const_model.h \
		multi_exp_Asqr_vec_BC_const_model.h \
		multi_exp_Asqr_expE_vec_BC_const_model.h \
		multi_alt_exp_vec_BC_const_model.h \
		multi_alt_exp_expE_vec_BC_const_model.h \
		multi_alt_exp_Asqr_vec_BC_const_model.h \
		multi_alt_exp_Asqr_expE_vec_BC_const_model.h \
		threept_multi_alt_exp_expE_model.h \
		threept_multi_alt_exp_expE_vec_model.h \
		threept_multi_alt_exp_model.h \
		threept_multi_alt_exp_vec_model.h \
		threept_multi_exp_expE_model.h \
		threept_multi_exp_expE_vec_model.h \
		threept_multi_exp_model.h \
		threept_multi_exp_vec_model.h \
		settingsmap.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o main.o main.cpp

settingsmap.o: settingsmap.cpp settingsmap.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o settingsmap.o settingsmap.cpp

binary_io.o: binary_io.cpp binary_io.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o binary_io.o binary_io.cpp

parse_model.o: parse_model.cpp parse_model.h \
		abstract_model.h \
		parser.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o parse_model.o parse_model.cpp

multi_alt_exp_Asqr_BC_model.o: multi_alt_exp_Asqr_BC_model.cpp multi_alt_exp_Asqr_BC_model.h \
		abstract_model.h \
		multi_alt_exp_Asqr_model.h \
		multi_exp_Asqr_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_Asqr_BC_model.o multi_alt_exp_Asqr_BC_model.cpp

multi_alt_exp_Asqr_expE_BC_model.o: multi_alt_exp_Asqr_expE_BC_model.cpp multi_alt_exp_Asqr_expE_BC_model.h \
		abstract_model.h \
		multi_alt_exp_Asqr_expE_model.h \
		multi_exp_Asqr_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_Asqr_expE_BC_model.o multi_alt_exp_Asqr_expE_BC_model.cpp

multi_alt_exp_Asqr_expE_model.o: multi_alt_exp_Asqr_expE_model.cpp multi_alt_exp_Asqr_expE_model.h \
		abstract_model.h \
		multi_exp_Asqr_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_Asqr_expE_model.o multi_alt_exp_Asqr_expE_model.cpp

multi_alt_exp_Asqr_expE_vec_BC_model.o: multi_alt_exp_Asqr_expE_vec_BC_model.cpp multi_alt_exp_Asqr_expE_vec_BC_model.h \
		abstract_model.h \
		multi_alt_exp_Asqr_expE_vec_model.h \
		multi_alt_exp_Asqr_expE_model.h \
		multi_exp_Asqr_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_Asqr_expE_vec_BC_model.o multi_alt_exp_Asqr_expE_vec_BC_model.cpp

multi_alt_exp_Asqr_expE_vec_model.o: multi_alt_exp_Asqr_expE_vec_model.cpp multi_alt_exp_Asqr_expE_vec_model.h \
		abstract_model.h \
		multi_alt_exp_Asqr_expE_model.h \
		multi_exp_Asqr_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_Asqr_expE_vec_model.o multi_alt_exp_Asqr_expE_vec_model.cpp

multi_alt_exp_Asqr_model.o: multi_alt_exp_Asqr_model.cpp multi_alt_exp_Asqr_model.h \
		abstract_model.h \
		multi_exp_Asqr_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_Asqr_model.o multi_alt_exp_Asqr_model.cpp

multi_alt_exp_Asqr_vec_BC_model.o: multi_alt_exp_Asqr_vec_BC_model.cpp multi_alt_exp_Asqr_vec_BC_model.h \
		abstract_model.h \
		multi_alt_exp_Asqr_vec_model.h \
		multi_alt_exp_Asqr_model.h \
		multi_exp_Asqr_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_Asqr_vec_BC_model.o multi_alt_exp_Asqr_vec_BC_model.cpp

multi_alt_exp_Asqr_vec_model.o: multi_alt_exp_Asqr_vec_model.cpp multi_alt_exp_Asqr_vec_model.h \
		abstract_model.h \
		multi_alt_exp_Asqr_model.h \
		multi_exp_Asqr_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_Asqr_vec_model.o multi_alt_exp_Asqr_vec_model.cpp

multi_alt_exp_BC_model.o: multi_alt_exp_BC_model.cpp multi_alt_exp_BC_model.h \
		abstract_model.h \
		multi_alt_exp_model.h \
		multi_exp_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_BC_model.o multi_alt_exp_BC_model.cpp

multi_alt_exp_expE_BC_model.o: multi_alt_exp_expE_BC_model.cpp multi_alt_exp_expE_BC_model.h \
		abstract_model.h \
		multi_alt_exp_expE_model.h \
		multi_exp_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_expE_BC_model.o multi_alt_exp_expE_BC_model.cpp

multi_alt_exp_expE_mat_model.o: multi_alt_exp_expE_mat_model.cpp multi_alt_exp_expE_mat_model.h \
		abstract_model.h \
		multi_exp_expE_mat_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_expE_mat_model.o multi_alt_exp_expE_mat_model.cpp

multi_alt_exp_expE_model.o: multi_alt_exp_expE_model.cpp multi_alt_exp_expE_model.h \
		abstract_model.h \
		multi_exp_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_expE_model.o multi_alt_exp_expE_model.cpp

multi_alt_exp_expE_nonsym_mat_model.o: multi_alt_exp_expE_nonsym_mat_model.cpp multi_alt_exp_expE_nonsym_mat_model.h \
		abstract_model.h \
		multi_exp_expE_nonsym_mat_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_expE_nonsym_mat_model.o multi_alt_exp_expE_nonsym_mat_model.cpp

multi_alt_exp_expE_vec_BC_model.o: multi_alt_exp_expE_vec_BC_model.cpp multi_alt_exp_expE_vec_BC_model.h \
		abstract_model.h \
		multi_alt_exp_expE_vec_model.h \
		multi_alt_exp_expE_model.h \
		multi_exp_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_expE_vec_BC_model.o multi_alt_exp_expE_vec_BC_model.cpp

multi_alt_exp_expE_vec_model.o: multi_alt_exp_expE_vec_model.cpp multi_alt_exp_expE_vec_model.h \
		abstract_model.h \
		multi_alt_exp_expE_model.h \
		multi_exp_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_expE_vec_model.o multi_alt_exp_expE_vec_model.cpp

multi_alt_exp_mat_model.o: multi_alt_exp_mat_model.cpp multi_alt_exp_mat_model.h \
		abstract_model.h \
		multi_exp_mat_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_mat_model.o multi_alt_exp_mat_model.cpp

multi_alt_exp_model.o: multi_alt_exp_model.cpp multi_alt_exp_model.h \
		abstract_model.h \
		multi_exp_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_model.o multi_alt_exp_model.cpp

multi_alt_exp_nonsym_mat_model.o: multi_alt_exp_nonsym_mat_model.cpp multi_alt_exp_nonsym_mat_model.h \
		abstract_model.h \
		multi_exp_nonsym_mat_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_nonsym_mat_model.o multi_alt_exp_nonsym_mat_model.cpp

multi_alt_exp_vec_BC_model.o: multi_alt_exp_vec_BC_model.cpp multi_alt_exp_vec_BC_model.h \
		abstract_model.h \
		multi_alt_exp_vec_model.h \
		multi_alt_exp_model.h \
		multi_exp_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_vec_BC_model.o multi_alt_exp_vec_BC_model.cpp

multi_alt_exp_vec_model.o: multi_alt_exp_vec_model.cpp multi_alt_exp_vec_model.h \
		abstract_model.h \
		multi_alt_exp_model.h \
		multi_exp_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_vec_model.o multi_alt_exp_vec_model.cpp

multi_exp_Asqr_BC_model.o: multi_exp_Asqr_BC_model.cpp multi_exp_Asqr_BC_model.h \
		abstract_model.h \
		multi_exp_Asqr_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_Asqr_BC_model.o multi_exp_Asqr_BC_model.cpp

multi_exp_Asqr_expE_BC_model.o: multi_exp_Asqr_expE_BC_model.cpp multi_exp_Asqr_expE_BC_model.h \
		abstract_model.h \
		multi_exp_Asqr_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_Asqr_expE_BC_model.o multi_exp_Asqr_expE_BC_model.cpp

multi_exp_Asqr_expE_model.o: multi_exp_Asqr_expE_model.cpp multi_exp_Asqr_expE_model.h \
		abstract_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_Asqr_expE_model.o multi_exp_Asqr_expE_model.cpp

multi_exp_Asqr_expE_vec_BC_model.o: multi_exp_Asqr_expE_vec_BC_model.cpp multi_exp_Asqr_expE_vec_BC_model.h \
		abstract_model.h \
		multi_exp_Asqr_expE_vec_model.h \
		multi_exp_Asqr_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_Asqr_expE_vec_BC_model.o multi_exp_Asqr_expE_vec_BC_model.cpp

multi_exp_Asqr_expE_vec_model.o: multi_exp_Asqr_expE_vec_model.cpp multi_exp_Asqr_expE_vec_model.h \
		abstract_model.h \
		multi_exp_Asqr_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_Asqr_expE_vec_model.o multi_exp_Asqr_expE_vec_model.cpp

multi_exp_Asqr_model.o: multi_exp_Asqr_model.cpp multi_exp_Asqr_model.h \
		abstract_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_Asqr_model.o multi_exp_Asqr_model.cpp

multi_exp_Asqr_vec_BC_model.o: multi_exp_Asqr_vec_BC_model.cpp multi_exp_Asqr_vec_BC_model.h \
		abstract_model.h \
		multi_exp_Asqr_vec_model.h \
		multi_exp_Asqr_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_Asqr_vec_BC_model.o multi_exp_Asqr_vec_BC_model.cpp

multi_exp_Asqr_vec_model.o: multi_exp_Asqr_vec_model.cpp multi_exp_Asqr_vec_model.h \
		abstract_model.h \
		multi_exp_Asqr_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_Asqr_vec_model.o multi_exp_Asqr_vec_model.cpp

multi_exp_BC_model.o: multi_exp_BC_model.cpp multi_exp_BC_model.h \
		abstract_model.h \
		multi_exp_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_BC_model.o multi_exp_BC_model.cpp

multi_exp_expE_BC_model.o: multi_exp_expE_BC_model.cpp multi_exp_expE_BC_model.h \
		abstract_model.h \
		multi_exp_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_expE_BC_model.o multi_exp_expE_BC_model.cpp

multi_exp_expE_mat_model.o: multi_exp_expE_mat_model.cpp multi_exp_expE_mat_model.h \
		abstract_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_expE_mat_model.o multi_exp_expE_mat_model.cpp

multi_exp_expE_mat_II_model.o: multi_exp_expE_mat_II_model.cpp multi_exp_expE_mat_II_model.h \
		abstract_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_expE_mat_II_model.o multi_exp_expE_mat_II_model.cpp

multi_exp_expE_mat_upper_model.o: multi_exp_expE_mat_upper_model.cpp multi_exp_expE_mat_upper_model.h \
		abstract_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_expE_mat_upper_model.o multi_exp_expE_mat_upper_model.cpp

multi_exp_expE_mat_II_upper_model.o: multi_exp_expE_mat_II_upper_model.cpp multi_exp_expE_mat_II_upper_model.h \
		abstract_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_expE_mat_II_upper_model.o multi_exp_expE_mat_II_upper_model.cpp

multi_exp_expE_model.o: multi_exp_expE_model.cpp multi_exp_expE_model.h \
		abstract_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_expE_model.o multi_exp_expE_model.cpp

multi_exp_expE_nonsym_mat_model.o: multi_exp_expE_nonsym_mat_model.cpp multi_exp_expE_nonsym_mat_model.h \
		abstract_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_expE_nonsym_mat_model.o multi_exp_expE_nonsym_mat_model.cpp

multi_exp_expE_vec_BC_model.o: multi_exp_expE_vec_BC_model.cpp multi_exp_expE_vec_BC_model.h \
		abstract_model.h \
		multi_exp_expE_vec_model.h \
		multi_exp_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_expE_vec_BC_model.o multi_exp_expE_vec_BC_model.cpp

multi_exp_expE_vec_model.o: multi_exp_expE_vec_model.cpp multi_exp_expE_vec_model.h \
		abstract_model.h \
		multi_exp_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_expE_vec_model.o multi_exp_expE_vec_model.cpp

multi_exp_mat_model.o: multi_exp_mat_model.cpp multi_exp_mat_model.h \
		abstract_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_mat_model.o multi_exp_mat_model.cpp

multi_exp_mat_II_model.o: multi_exp_mat_II_model.cpp multi_exp_mat_II_model.h \
		abstract_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_mat_II_model.o multi_exp_mat_II_model.cpp

multi_exp_mat_upper_model.o: multi_exp_mat_upper_model.cpp multi_exp_mat_upper_model.h \
		abstract_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_mat_upper_model.o multi_exp_mat_upper_model.cpp

multi_exp_mat_II_upper_model.o: multi_exp_mat_II_upper_model.cpp multi_exp_mat_II_upper_model.h \
		abstract_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_mat_II_upper_model.o multi_exp_mat_II_upper_model.cpp

multi_exp_model.o: multi_exp_model.cpp multi_exp_model.h \
		abstract_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_model.o multi_exp_model.cpp

multi_exp_nonsym_mat_model.o: multi_exp_nonsym_mat_model.cpp multi_exp_nonsym_mat_model.h \
		abstract_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_nonsym_mat_model.o multi_exp_nonsym_mat_model.cpp

multi_exp_vec_BC_model.o: multi_exp_vec_BC_model.cpp multi_exp_vec_BC_model.h \
		abstract_model.h \
		multi_exp_vec_model.h \
		multi_exp_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_vec_BC_model.o multi_exp_vec_BC_model.cpp

multi_exp_vec_model.o: multi_exp_vec_model.cpp multi_exp_vec_model.h \
		abstract_model.h \
		multi_exp_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_vec_model.o multi_exp_vec_model.cpp

multi_exp_const_model.o: multi_exp_const_model.cpp multi_exp_const_model.h \
		abstract_model.h \
		multi_exp_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_const_model.o multi_exp_const_model.cpp

multi_exp_expE_const_model.o: multi_exp_expE_const_model.cpp multi_exp_expE_const_model.h \
		abstract_model.h \
		multi_exp_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_expE_const_model.o multi_exp_expE_const_model.cpp

multi_exp_Asqr_const_model.o: multi_exp_Asqr_const_model.cpp multi_exp_Asqr_const_model.h \
		abstract_model.h \
		multi_exp_Asqr_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_Asqr_const_model.o multi_exp_Asqr_const_model.cpp

multi_exp_Asqr_expE_const_model.o: multi_exp_Asqr_expE_const_model.cpp multi_exp_Asqr_expE_const_model.h \
		abstract_model.h \
		multi_exp_Asqr_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_Asqr_expE_const_model.o multi_exp_Asqr_expE_const_model.cpp

multi_alt_exp_const_model.o: multi_alt_exp_const_model.cpp multi_alt_exp_const_model.h \
		abstract_model.h \
		multi_alt_exp_model.h \
		multi_exp_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_const_model.o multi_alt_exp_const_model.cpp

multi_alt_exp_expE_const_model.o: multi_alt_exp_expE_const_model.cpp multi_alt_exp_expE_const_model.h \
		abstract_model.h \
		multi_alt_exp_expE_model.h \
		multi_exp_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_expE_const_model.o multi_alt_exp_expE_const_model.cpp

multi_alt_exp_Asqr_const_model.o: multi_alt_exp_Asqr_const_model.cpp multi_alt_exp_Asqr_const_model.h \
		abstract_model.h \
		multi_alt_exp_Asqr_model.h \
		multi_exp_Asqr_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_Asqr_const_model.o multi_alt_exp_Asqr_const_model.cpp

multi_alt_exp_Asqr_expE_const_model.o: multi_alt_exp_Asqr_expE_const_model.cpp multi_alt_exp_Asqr_expE_const_model.h \
		abstract_model.h \
		multi_alt_exp_Asqr_expE_model.h \
		multi_exp_Asqr_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_Asqr_expE_const_model.o multi_alt_exp_Asqr_expE_const_model.cpp

multi_exp_vec_const_model.o: multi_exp_vec_const_model.cpp multi_exp_vec_const_model.h \
		abstract_model.h \
		multi_exp_vec_model.h \
		multi_exp_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_vec_const_model.o multi_exp_vec_const_model.cpp

multi_exp_expE_vec_const_model.o: multi_exp_expE_vec_const_model.cpp multi_exp_expE_vec_const_model.h \
		abstract_model.h \
		multi_exp_expE_vec_model.h \
		multi_exp_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_expE_vec_const_model.o multi_exp_expE_vec_const_model.cpp

multi_exp_Asqr_vec_const_model.o: multi_exp_Asqr_vec_const_model.cpp multi_exp_Asqr_vec_const_model.h \
		abstract_model.h \
		multi_exp_Asqr_vec_model.h \
		multi_exp_Asqr_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_Asqr_vec_const_model.o multi_exp_Asqr_vec_const_model.cpp

multi_exp_Asqr_expE_vec_const_model.o: multi_exp_Asqr_expE_vec_const_model.cpp multi_exp_Asqr_expE_vec_const_model.h \
		abstract_model.h \
		multi_exp_Asqr_expE_vec_model.h \
		multi_exp_Asqr_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_Asqr_expE_vec_const_model.o multi_exp_Asqr_expE_vec_const_model.cpp

multi_alt_exp_vec_const_model.o: multi_alt_exp_vec_const_model.cpp multi_alt_exp_vec_const_model.h \
		abstract_model.h \
		multi_alt_exp_vec_model.h \
		multi_alt_exp_model.h \
		multi_exp_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_vec_const_model.o multi_alt_exp_vec_const_model.cpp

multi_alt_exp_expE_vec_const_model.o: multi_alt_exp_expE_vec_const_model.cpp multi_alt_exp_expE_vec_const_model.h \
		abstract_model.h \
		multi_alt_exp_expE_vec_model.h \
		multi_alt_exp_expE_model.h \
		multi_exp_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_expE_vec_const_model.o multi_alt_exp_expE_vec_const_model.cpp

multi_alt_exp_Asqr_vec_const_model.o: multi_alt_exp_Asqr_vec_const_model.cpp multi_alt_exp_Asqr_vec_const_model.h \
		abstract_model.h \
		multi_alt_exp_Asqr_vec_model.h \
		multi_alt_exp_Asqr_model.h \
		multi_exp_Asqr_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_Asqr_vec_const_model.o multi_alt_exp_Asqr_vec_const_model.cpp

multi_alt_exp_Asqr_expE_vec_const_model.o: multi_alt_exp_Asqr_expE_vec_const_model.cpp multi_alt_exp_Asqr_expE_vec_const_model.h \
		abstract_model.h \
		multi_alt_exp_Asqr_expE_vec_model.h \
		multi_alt_exp_Asqr_expE_model.h \
		multi_exp_Asqr_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_Asqr_expE_vec_const_model.o multi_alt_exp_Asqr_expE_vec_const_model.cpp

multi_exp_BC_const_model.o: multi_exp_BC_const_model.cpp multi_exp_BC_const_model.h \
		abstract_model.h \
		multi_exp_BC_model.h \
		multi_exp_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_BC_const_model.o multi_exp_BC_const_model.cpp

multi_exp_expE_BC_const_model.o: multi_exp_expE_BC_const_model.cpp multi_exp_expE_BC_const_model.h \
		abstract_model.h \
		multi_exp_expE_BC_model.h \
		multi_exp_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_expE_BC_const_model.o multi_exp_expE_BC_const_model.cpp

multi_exp_Asqr_BC_const_model.o: multi_exp_Asqr_BC_const_model.cpp multi_exp_Asqr_BC_const_model.h \
		abstract_model.h \
		multi_exp_Asqr_BC_model.h \
		multi_exp_Asqr_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_Asqr_BC_const_model.o multi_exp_Asqr_BC_const_model.cpp

multi_exp_Asqr_expE_BC_const_model.o: multi_exp_Asqr_expE_BC_const_model.cpp multi_exp_Asqr_expE_BC_const_model.h \
		abstract_model.h \
		multi_exp_Asqr_expE_BC_model.h \
		multi_exp_Asqr_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_Asqr_expE_BC_const_model.o multi_exp_Asqr_expE_BC_const_model.cpp

multi_alt_exp_BC_const_model.o: multi_alt_exp_BC_const_model.cpp multi_alt_exp_BC_const_model.h \
		abstract_model.h \
		multi_alt_exp_BC_model.h \
		multi_alt_exp_model.h \
		multi_exp_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_BC_const_model.o multi_alt_exp_BC_const_model.cpp

multi_alt_exp_expE_BC_const_model.o: multi_alt_exp_expE_BC_const_model.cpp multi_alt_exp_expE_BC_const_model.h \
		abstract_model.h \
		multi_alt_exp_expE_BC_model.h \
		multi_alt_exp_expE_model.h \
		multi_exp_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_expE_BC_const_model.o multi_alt_exp_expE_BC_const_model.cpp

multi_alt_exp_Asqr_BC_const_model.o: multi_alt_exp_Asqr_BC_const_model.cpp multi_alt_exp_Asqr_BC_const_model.h \
		abstract_model.h \
		multi_alt_exp_Asqr_BC_model.h \
		multi_alt_exp_Asqr_model.h \
		multi_exp_Asqr_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_Asqr_BC_const_model.o multi_alt_exp_Asqr_BC_const_model.cpp

multi_alt_exp_Asqr_expE_BC_const_model.o: multi_alt_exp_Asqr_expE_BC_const_model.cpp multi_alt_exp_Asqr_expE_BC_const_model.h \
		abstract_model.h \
		multi_alt_exp_Asqr_expE_BC_model.h \
		multi_alt_exp_Asqr_expE_model.h \
		multi_exp_Asqr_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_Asqr_expE_BC_const_model.o multi_alt_exp_Asqr_expE_BC_const_model.cpp

multi_exp_vec_BC_const_model.o: multi_exp_vec_BC_const_model.cpp multi_exp_vec_BC_const_model.h \
		abstract_model.h \
		multi_exp_vec_BC_model.h \
		multi_exp_vec_model.h \
		multi_exp_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_vec_BC_const_model.o multi_exp_vec_BC_const_model.cpp

multi_exp_expE_vec_BC_const_model.o: multi_exp_expE_vec_BC_const_model.cpp multi_exp_expE_vec_BC_const_model.h \
		abstract_model.h \
		multi_exp_expE_vec_BC_model.h \
		multi_exp_expE_vec_model.h \
		multi_exp_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_expE_vec_BC_const_model.o multi_exp_expE_vec_BC_const_model.cpp

multi_exp_Asqr_vec_BC_const_model.o: multi_exp_Asqr_vec_BC_const_model.cpp multi_exp_Asqr_vec_BC_const_model.h \
		abstract_model.h \
		multi_exp_Asqr_vec_BC_model.h \
		multi_exp_Asqr_vec_model.h \
		multi_exp_Asqr_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_Asqr_vec_BC_const_model.o multi_exp_Asqr_vec_BC_const_model.cpp

multi_exp_Asqr_expE_vec_BC_const_model.o: multi_exp_Asqr_expE_vec_BC_const_model.cpp multi_exp_Asqr_expE_vec_BC_const_model.h \
		abstract_model.h \
		multi_exp_Asqr_expE_vec_BC_model.h \
		multi_exp_Asqr_expE_vec_model.h \
		multi_exp_Asqr_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_exp_Asqr_expE_vec_BC_const_model.o multi_exp_Asqr_expE_vec_BC_const_model.cpp

multi_alt_exp_vec_BC_const_model.o: multi_alt_exp_vec_BC_const_model.cpp multi_alt_exp_vec_BC_const_model.h \
		abstract_model.h \
		multi_alt_exp_vec_BC_model.h \
		multi_alt_exp_vec_model.h \
		multi_alt_exp_model.h \
		multi_exp_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_vec_BC_const_model.o multi_alt_exp_vec_BC_const_model.cpp

multi_alt_exp_expE_vec_BC_const_model.o: multi_alt_exp_expE_vec_BC_const_model.cpp multi_alt_exp_expE_vec_BC_const_model.h \
		abstract_model.h \
		multi_alt_exp_expE_vec_BC_model.h \
		multi_alt_exp_expE_vec_model.h \
		multi_alt_exp_expE_model.h \
		multi_exp_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_expE_vec_BC_const_model.o multi_alt_exp_expE_vec_BC_const_model.cpp

multi_alt_exp_Asqr_vec_BC_const_model.o: multi_alt_exp_Asqr_vec_BC_const_model.cpp multi_alt_exp_Asqr_vec_BC_const_model.h \
		abstract_model.h \
		multi_alt_exp_Asqr_vec_BC_model.h \
		multi_alt_exp_Asqr_vec_model.h \
		multi_alt_exp_Asqr_model.h \
		multi_exp_Asqr_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_Asqr_vec_BC_const_model.o multi_alt_exp_Asqr_vec_BC_const_model.cpp

multi_alt_exp_Asqr_expE_vec_BC_const_model.o: multi_alt_exp_Asqr_expE_vec_BC_const_model.cpp multi_alt_exp_Asqr_expE_vec_BC_const_model.h \
		abstract_model.h \
		multi_alt_exp_Asqr_expE_vec_BC_model.h \
		multi_alt_exp_Asqr_expE_vec_model.h \
		multi_alt_exp_Asqr_expE_model.h \
		multi_exp_Asqr_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multi_alt_exp_Asqr_expE_vec_BC_const_model.o multi_alt_exp_Asqr_expE_vec_BC_const_model.cpp

threept_multi_alt_exp_expE_model.o: threept_multi_alt_exp_expE_model.cpp threept_multi_alt_exp_expE_model.h \
		abstract_model.h \
		threept_multi_exp_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o threept_multi_alt_exp_expE_model.o threept_multi_alt_exp_expE_model.cpp

threept_multi_alt_exp_expE_vec_model.o: threept_multi_alt_exp_expE_vec_model.cpp threept_multi_alt_exp_expE_vec_model.h \
		abstract_model.h \
		threept_multi_alt_exp_expE_model.h \
		threept_multi_exp_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o threept_multi_alt_exp_expE_vec_model.o threept_multi_alt_exp_expE_vec_model.cpp

threept_multi_alt_exp_model.o: threept_multi_alt_exp_model.cpp threept_multi_alt_exp_model.h \
		abstract_model.h \
		threept_multi_exp_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o threept_multi_alt_exp_model.o threept_multi_alt_exp_model.cpp

threept_multi_alt_exp_vec_model.o: threept_multi_alt_exp_vec_model.cpp threept_multi_alt_exp_vec_model.h \
		abstract_model.h \
		threept_multi_alt_exp_model.h \
		threept_multi_exp_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o threept_multi_alt_exp_vec_model.o threept_multi_alt_exp_vec_model.cpp

threept_multi_exp_expE_model.o: threept_multi_exp_expE_model.cpp threept_multi_exp_expE_model.h \
		abstract_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o threept_multi_exp_expE_model.o threept_multi_exp_expE_model.cpp

threept_multi_exp_expE_vec_model.o: threept_multi_exp_expE_vec_model.cpp threept_multi_exp_expE_vec_model.h \
		abstract_model.h \
		threept_multi_exp_expE_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o threept_multi_exp_expE_vec_model.o threept_multi_exp_expE_vec_model.cpp

threept_multi_exp_model.o: threept_multi_exp_model.cpp threept_multi_exp_model.h \
		abstract_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o threept_multi_exp_model.o threept_multi_exp_model.cpp

threept_multi_exp_vec_model.o: threept_multi_exp_vec_model.cpp threept_multi_exp_vec_model.h \
		abstract_model.h \
		threept_multi_exp_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o threept_multi_exp_vec_model.o threept_multi_exp_vec_model.cpp

gaussian_prior.o: gaussian_prior.cpp gaussian_prior.h \
		abstract_prior.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o gaussian_prior.o gaussian_prior.cpp

fitter.o: fitter.cpp fitter.h \
		abstract_model.h \
		abstract_prior.h \
		chisqr_extra_term.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o fitter.o fitter.cpp

parser.o: parser.cpp parser.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o parser.o parser.cpp

chisqr_extra_term.o: chisqr_extra_term.cpp chisqr_extra_term.h \
		parser.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o chisqr_extra_term.o chisqr_extra_term.cpp

.PHONY : clean
clean :
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) *~
