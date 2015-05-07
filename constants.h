#define INCREMENT   10

#define VARS            4       //number of variables per particle (3 states and a W)
#define NUM_PARTICLES   10000    //number of particles in our distribution
#define NUM_PARTICLES_F 10000.0   //float version of particles for average

#define X            0       //labels for variable positions in particle array
#define Y            1
#define Q            2
#define W            3

#define PI              3.14159

#define MAX_X           300.0   //max values for position and heading given the environment
#define MAX_Y           200.0
#define MAX_Q           180.0

#define MAX_X_INT           300
#define MAX_Y_INT           200
#define MAX_Q_INT           180

#define ENCODER_NOISE   0.5    //amount of noise in encoders (as a fraction of the 
#define POSITION_NOISE  10.0   //amount of noise in position in mm

#define MAX_ANGLE 90
#define ANGLE_STEP 15
#define OBS_SIZE_X 40
#define OBS_SIZE_Y 40
#define SONAR_SIGMA 2
#define MAX_RANGE 253
#define CUTOFF_RANGE 100

#define FIRST_STEP 10
#define SECOND_STEP 5

#define NUM_OBSTACLES 7

#define CONVERGED 10
