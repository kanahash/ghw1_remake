NAME = ghw1

CXX = g++
CXXFLAGS = -I./include -mcmodel=large -march=native -O2 -fopenmp
# LDFLAGSに -fopenmp を追加します
LDFLAGS = -fopenmp -lm -lfftw3_omp -lfftw3 
LIBS_PATH = /usr/lib64

SRCS = srcs/main.cpp\
        srcs/init/init.cpp\
        srcs/init/init_output_and_arrays.cpp\
        srcs/init/init_parameters.cpp\
        srcs/init/init_add_turb.cpp\
        srcs/setting/set_init_density.cpp\
        srcs/setting/set_time.cpp\
        srcs/setting/set_trace.cpp\
        srcs/calculation/init_potential.cpp\
        srcs/calculation/density_time_evolution_step.cpp\
        srcs/calculation/gyro_average_calculate.cpp\
        srcs/calculation/solve_polarization_equation.cpp\
        srcs/calculation/gyro_shielded_potential_calculate.cpp\
        srcs/calculation/precompute_fftw_coefficients.cpp\
        srcs/calculation/FFT.cpp\
        srcs/calculation/energy_and_transport_calculate.cpp\
        srcs/calculation/linear_growth_and_frequency_calculate.cpp\
        srcs/calculation/calculate_and_write_ky_spectrum.cpp\
        srcs/calculation/calculate_and_write_energy_ky_spectrum.cpp\
        srcs/calculation/calculate_and_write_kx_spectrum.cpp\
        srcs/calculation/diagnose.cpp\
        srcs/calculation/arakawa.cpp\
        srcs/calculation/arakaw4.cpp\
        srcs/calculation/laplace.cpp\
        srcs/calculation/poisson.cpp\
        srcs/utils/utils_time.cpp\
        srcs/utils/utils_history.cpp\
        srcs/utils/utils_data.cpp

OBJS = $(SRCS:.cpp=.o)

all: $(NAME)

$(NAME): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS) -L$(LIBS_PATH)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJS)

fclean: clean
	$(RM) $(NAME)

re: fclean all

.PHONY: all clean fclean re
