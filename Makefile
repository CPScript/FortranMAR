# Compiler selection. can be overridden with make FC=ifort
FC ?= gfortran
RSCRIPT ?= Rscript

FFLAGS_COMMON = -std=f2008 -pedantic -Wall
FFLAGS_DEBUG = $(FFLAGS_COMMON) -g -fbacktrace -fbounds-check -O0
FFLAGS_RELEASE = $(FFLAGS_COMMON) -O3 -march=native -funroll-loops

ifeq ($(FC), ifort)
    FFLAGS_DEBUG = -stand f08 -warn all -g -traceback -check bounds -O0
    FFLAGS_RELEASE = -stand f08 -O3 -xHost -unroll
endif

FFLAGS ?= $(FFLAGS_RELEASE)

SRC = projectile_motion.f90
EXE = projectile_motion
R_SCRIPT = trajectory_analysis.R

OBJ = $(SRC:.f90=.o)

.PHONY: all
all: $(EXE)

$(EXE): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

.PHONY: debug
debug: FFLAGS = $(FFLAGS_DEBUG)
debug: $(EXE)

.PHONY: release
release: FFLAGS = $(FFLAGS_RELEASE)
release: $(EXE)

.PHONY: run
run: $(EXE)
	./$(EXE)

.PHONY: analyze
analyze: $(EXE)
	@echo "Running Fortran simulation..."
	./$(EXE)
	@echo "Running R analysis..."
	$(RSCRIPT) $(R_SCRIPT)
	@echo "Analysis complete. Check generated files."

.PHONY: profile
profile: FFLAGS = $(FFLAGS_RELEASE) -pg
profile: $(EXE)
	./$(EXE)
	gprof $(EXE) gmon.out > profile_report.txt
	@echo "Profile report generated: profile_report.txt"

.PHONY: memcheck
memcheck: debug
	valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all ./$(EXE)

.PHONY: clean
clean:
	rm -f $(OBJ) $(EXE) *.mod gmon.out profile_report.txt

.PHONY: distclean
distclean: clean
	rm -f trajectory_data.csv processed_trajectory_data.csv
	rm -f trajectory_analysis.png trajectory_comparison.png
	rm -f interactive_trajectory.html

.PHONY: install-deps
install-deps:
	$(RSCRIPT) -e "install.packages(c('ggplot2', 'dplyr', 'gridExtra', 'viridis', 'plotly', 'tidyr', 'htmlwidgets'), repos='https://cran.r-project.org')"

.PHONY: check-deps
check-deps:
	$(RSCRIPT) -e "required <- c('ggplot2', 'dplyr', 'gridExtra', 'viridis', 'plotly', 'tidyr'); installed <- installed.packages()[,'Package']; missing <- required[!required %in% installed]; if(length(missing) > 0) { cat('Missing packages:', paste(missing, collapse=', '), '\n'); quit(status=1) } else { cat('All required packages installed.\n') }"

.PHONY: benchmark
benchmark: release
	@echo "Running performance benchmark..."
	@time ./$(EXE) > /dev/null
	@echo "Benchmark complete."

.PHONY: docs
docs:
	@echo "=== Project Documentation ===" > README.txt
	@echo "Projectile Motion Simulation with Air Resistance" >> README.txt
	@echo "" >> README.txt
	@echo "Build commands:" >> README.txt
	@echo "  make          - Build release version" >> README.txt
	@echo "  make debug    - Build debug version" >> README.txt
	@echo "  make run      - Run simulation" >> README.txt
	@echo "  make analyze  - Run full analysis pipeline" >> README.txt
	@echo "  make clean    - Clean build files" >> README.txt
	@echo "" >> README.txt
	@echo "Files generated:" >> README.txt
	@echo "  trajectory_data.csv - Raw simulation data" >> README.txt
	@echo "  *.png - Analysis plots" >> README.txt
	@echo "  *.html - Interactive visualizations" >> README.txt
	@echo "" >> README.txt
	@echo "Generated: $$(date)" >> README.txt
	@echo "Documentation created: README.txt"

.PHONY: help
help:
	@echo "Available targets:"
	@echo "  all         - Build executable (default)"
	@echo "  debug       - Build with debugging symbols"
	@echo "  release     - Build optimized version"
	@echo "  run         - Execute simulation"
	@echo "  analyze     - Run simulation and R analysis"
	@echo "  profile     - Build and profile executable"
	@echo "  memcheck    - Run with valgrind memory checker"
	@echo "  benchmark   - Time execution performance"
	@echo "  install-deps- Install required R packages"
	@echo "  check-deps  - Check R package dependencies"
	@echo "  docs        - Generate documentation"
	@echo "  clean       - Remove build artifacts"
	@echo "  distclean   - Remove all generated files"
	@echo "  help        - Show this help message"
	@echo ""
	@echo "Variables:"
	@echo "  FC          - Fortran compiler (default: gfortran)"
	@echo "  RSCRIPT     - R script executor (default: Rscript)"
	@echo "  FFLAGS      - Compiler flags"
