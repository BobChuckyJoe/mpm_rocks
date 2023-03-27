// Simulation parameters
const double BOX_SIZE = 10.0;
const long GRID_LENGTH = 4;
const float GRID_SPACING = BOX_SIZE / GRID_LENGTH;
const double DELTA_T = 0.001;
const long N_PARTICLES = 4;
const long N_ITERATIONS = 10;

// Consituency model parameters
const double CRITICAL_COMPRESSION = 2.5e-2;
const double CRITICAL_STRETCH = 7.5e-3;
const double HARDENING_COEFFICIENT = 10.0;
const double INITIAL_DENSITY = 4e2;
const double YOUNGS_MODULUS = 1.4e5;
const double POISSONS_RATIO = 0.2;

const double MEW_NOUGHT = YOUNGS_MODULUS / (2 * (1 + POISSONS_RATIO));
const double LAMBDA_NOUGHT = YOUNGS_MODULUS * POISSONS_RATIO / ( ( 1 + POISSONS_RATIO ) * (1.0 - 2.0 * POISSONS_RATIO));