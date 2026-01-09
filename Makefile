CXX      := g++
CXXFLAGS := -std=c++20 -Wall -Wextra -O2

TARGETS := test_v2d test_particle main

# Default: build everything
all: $(TARGETS)

# ------------------------------------------------
# Build V2D test
# ------------------------------------------------
test_v2d: test_v2d.cpp v2d.h
	$(CXX) $(CXXFLAGS) test_v2d.cpp -o $@

# ------------------------------------------------
# Build Particle test
# ------------------------------------------------
test_particle: test_particle.cpp particle.h v2d.h
	$(CXX) $(CXXFLAGS) test_particle.cpp -o $@

# ------------------------------------------------
# Build main executable
# ------------------------------------------------
main: main.cpp particle.h v2d.h
	$(CXX) $(CXXFLAGS) main.cpp -o $@

# ------------------------------------------------
# Run tests only
# ------------------------------------------------
test: test_v2d test_particle
	@echo "Running V2D tests..."
	./test_v2d
	@echo "Running Particle tests..."
	./test_particle

# ------------------------------------------------
# Clean
# ------------------------------------------------
clean:
	rm -f $(TARGETS)
