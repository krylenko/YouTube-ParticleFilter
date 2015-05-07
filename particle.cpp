/************************** Particle Filter Code *****************************/
/****************** Dan Ford, Gavin Ransom, Ramaprasad Kulkarni **************/
/****************** ECE 516 Robotics Fall 2009 *******************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <time.h>
#include "constants.h"          //header file with relevant constants

float x,y,q;                    //estimated position variables
float dp,dq;                  //movement inputs from robot simulator code
float sx,sy,sq;                 //variances of position variables
int converged;
float x_rand,y_rand,q_rand;     //random variables corrected to fit in environment
float rand_value;
//int w;                          //particle weight variable
int neighborhood;

float tx, ty, tq;
int goal_x,goal_y;

//array of particles with positions, headings, and Ws
float particle_map[VARS][NUM_PARTICLES];
int obstacles[2][NUM_OBSTACLES];

//array of sonar distance readings
int sonar[MAX_ANGLE*2/ANGLE_STEP+1];

//function prototypes
void particle_filter(void);
void init_dist(void);
void resample(void);
void print_file(int num);
void print_obst(void);
float dist(float x,float y);
float random_num(void);
float gaussian(void);
void update_position();
void weight_particles(void);
void print_init(void);
void generate_outputs(int num);
int obst_detect(float x,float y);
void hist(void);
void generateSonar();
void planPath();
void moveRobot();
void print_path(int num);
int sign(int v);
int getSonar(float x, float y, float q);
int getSonar2(float x, float y, float q);
float distance(float x0,float y0, float x1,float y1);

int main()
{
    //locate obstacles
	obstacles[X][0] = 50;      
	obstacles[Y][0] = 125;		
	obstacles[X][1] = 75;		
	obstacles[Y][1] = 50;		
	obstacles[X][2] = 125;		
	obstacles[Y][2] = 150;
	obstacles[X][3] = 150;
	obstacles[Y][3] = 25;
	obstacles[X][4] = 180;
	obstacles[Y][4] = 90;
	obstacles[X][5] = 250;
	obstacles[Y][5] = 125;
	obstacles[X][6] = 240;
	obstacles[Y][6] = 40;

    
    srand(time(NULL));    //seed random number generator
	print_obst();         //print obstacles to file

	//define goal location
    goal_x=280;
    goal_y=80;
    dp = dq = 0;
	sx = sy = sq = 100;
	converged = 0;
	FILE *out = fopen("output/goal.txt","w");
	fprintf(out,"%d\t%d",goal_x,goal_y);
	fclose(out);

    tx = 5.0f;
	ty = 80.0f;
	tq = 0.0f;	
	/*
	tx = random_num()*MAX_X;
	ty = random_num()*MAX_Y;
	tq = random_num()*MAX_Q*2;
	while(obst_detect(tx,ty))
	{
		tx = random_num()*MAX_X;
		ty = random_num()*MAX_Y;
	}*/

    //generate initial uniform distribution and save it to a file
	init_dist();
	print_init();

	FILE *outp = fopen("output/particle.txt","w");
	fprintf(outp,"%.2f\t%.2f\n", tx, ty);
	int i = 0;
	//move robot and check filter
	for(int j = 0; j < 200; ++j)
	//while( (abs(tx-goal_x)>1.0) || (abs(ty-goal_y)>1.0) )
	{
		generateSonar();
		particle_filter();
		generate_outputs(i);
        print_path(i++);
		fprintf(outp,"%.2f\t%.2f\n", tx, ty);
		planPath();
		moveRobot();
	}
	
	fclose(outp);  
    
	return 0;
}

//generate simulated noisy sonar readings
void generateSonar()
{
	for(int i = 0; i < MAX_ANGLE*2/ANGLE_STEP+1; ++i)
	{
		sonar[i] = (getSonar(tx,ty,i*ANGLE_STEP-MAX_ANGLE+tq) + gaussian()*SONAR_SIGMA);
		if( sonar[i]>CUTOFF_RANGE )
		{
			sonar[i] = CUTOFF_RANGE;
		}
	}
}

//move robot position
void moveRobot()
{
	tq += dq;
	if(tq > MAX_Q_INT)
   	{
   		tq -= 360;
   	}
	if(tq < -MAX_Q_INT)
	{
		tq += 360;
	}
	tx += dp*cos(tq*PI/180);
	ty += dp*sin(tq*PI/180);
}

//plan motion through map based on obstacles detected
void planPath()
{
	int dist = 2;
	dp = dist;
	if(!converged)
	{
		if(sonar[MAX_ANGLE/ANGLE_STEP] < 5)
		{
			dq = 90;
		}
		else
		{
			int edge[MAX_ANGLE/ANGLE_STEP*2+1];
			int count = 0;
			if(sonar[0] < CUTOFF_RANGE/2)
			{
				edge[count++] = -90;
			}
			for(int i = 1; i < MAX_ANGLE/ANGLE_STEP*2+1; ++i)
			{
				if( abs(sonar[i] - sonar[i-1]) > 20 )
				{
					edge[count++] = i*ANGLE_STEP - MAX_ANGLE;
				}
			}
			if(count < 3)
			{
				dq = 0;
			}
			else
			{
				dq = (edge[2] + edge[1])/2;
			}
		}
	}
	//use Tangent Bug algorithm to circumvent obstacles
	else
	{	
		dq = atan2(goal_y-y,goal_x-x)*180/PI - q;
		int dq_int = dq/ANGLE_STEP;
		int start = MAX_ANGLE/ANGLE_STEP + dq_int;
		if(sonar[start] < CUTOFF_RANGE)
		{
			int i = 1;
			int range;
			float left,right;
			float right_q,left_q;
			int current;

			while((i+dq_int)*ANGLE_STEP < MAX_ANGLE)
			{
				current = MAX_ANGLE/ANGLE_STEP+dq_int-i;
				if( (sonar[current] - sonar[start]) > 15 )
				{
					break;
				}
				++i;
			}
			range = sonar[MAX_ANGLE/ANGLE_STEP-(i-1)];
			right_q = q - (i-1+dq_int)*ANGLE_STEP;
			right = range + distance(goal_x,goal_y,x+range*cos(right_q*PI/180),y+range*sin(right_q*PI/180));
			
			i = 1;
			while((i+dq_int)*ANGLE_STEP < MAX_ANGLE)
			{
				current = MAX_ANGLE/ANGLE_STEP+dq_int+i;
				if( (sonar[current] - sonar[start]) > 15  )
				{
					break;
				}
				++i;
			}
			range = sonar[MAX_ANGLE/ANGLE_STEP+(i-1)];
			left_q = q + (i-1+dq_int)*ANGLE_STEP;
			left = range + distance(goal_x,goal_y,x+range*cos(left_q*PI/180),y+range*sin(left_q*PI/180));
			if(right < left)
			{
				dq = right_q;
			}
			else
			{
				dq = left_q;
			}
		}
	}
}

//main particle filtering function
void particle_filter(void)
{
	update_position();
	weight_particles();
 	resample();
}

//generate initial uniform distribution
void init_dist(void){
     
     int i;
     
     //set all particle weights identically
     for(i=0;i<NUM_PARTICLES;i++)
	 {
	     particle_map[W][i]=1/NUM_PARTICLES_F;
	     
	     //gen X          
		 particle_map[X][i]=random_num()*MAX_X;

	     //gen Y
	     particle_map[Y][i]=random_num()*MAX_Y;

	     //gen Q
		 particle_map[Q][i]= random_num()*MAX_Q*2;
	     
	     //check for obstacles
	     if(obst_detect(particle_map[X][i],particle_map[Y][i]) )
		 {
	     	--i;
			continue;
		 }
     }       
}

void resample(void){

        int j;
        float new_x[NUM_PARTICLES];         //init new x position array   
        float new_y[NUM_PARTICLES];         //init new y position array
        float new_q[NUM_PARTICLES];         //init new heading array
        float cdf[NUM_PARTICLES];           //init CDF array
        float ran[NUM_PARTICLES];           //init random number array

        //initialize CDF
        cdf[0] = particle_map[W][0];        
        for(int i = 1; i < NUM_PARTICLES; ++i)
        {  
                cdf[i] = cdf[i-1] + particle_map[W][i];
        }

        //initialize random array 
        for(int i = 0; i < NUM_PARTICLES; ++i)
        {
                ran[i] = rand()*cdf[NUM_PARTICLES-1]/RAND_MAX;
        }
		
		//make copies as many times as numbers occur in CDF
        for(int i = 0; i < NUM_PARTICLES; ++i)
        {
                j = 0;
                while(cdf[j] < ran[i])
                {
                        ++j;
                }
                new_x[i] = particle_map[X][j];
                new_y[i] = particle_map[Y][j];
                new_q[i] = particle_map[Q][j];
        }

        for(int i = 0; i < NUM_PARTICLES; ++i)
        {
                particle_map[X][i] = new_x[i];
                particle_map[Y][i] = new_y[i];
                particle_map[Q][i] = new_q[i];
                particle_map[W][i] = 1.0/NUM_PARTICLES_F;
        }

       
}

//generate random number scaled between 0 and 1
float random_num(void){
            
      return (rand()/(float)RAND_MAX);
      
}
 
//generate gaussian number for measurement noise     
float gaussian(void){
     float sum = 0.0;
	 int it = 10;
	 for (int i = 0; i < it; ++i)
	 {
		 sum += random_num()-.5;
	 }
     //return sqrt(-2.0 * log( (float)random_num()) ) * cos(2.0 * PI * random_num());
     return sum*1.0871;
}

//update positions of particles based on movement model that includes noise
//if particles leave the map, or a heading is less than 0 or greater than 2pi,
//generate a new random particle
void update_position()
{
    int i;
    for(i=0;i<NUM_PARTICLES;i++)
    {
	   particle_map[Q][i] += dq*(1 + gaussian()*ENCODER_NOISE);
	   if(particle_map[Q][i] > MAX_Q_INT)
	   {
	   		particle_map[Q][i] -= 360;
	   }
	   if(particle_map[Q][i] < -MAX_Q_INT)
	   {
	   		particle_map[Q][i] += 360;
	   }
	   particle_map[X][i] += dp*cos(particle_map[Q][i]*PI/180)*(1 + gaussian()*ENCODER_NOISE);  
	   particle_map[Y][i] += dp*sin(particle_map[Q][i]*PI/180)*(1 + gaussian()*ENCODER_NOISE);
    }
}

void weight_particles(void){
     //weight particles based on sonar matches          
     int i=0,j=0,diff=0;
     double discrep=0;
     int s_part[MAX_ANGLE*2/ANGLE_STEP+1];
    
    generateSonar();
    
    for( i=0;i<NUM_PARTICLES;i++)
	{
		if( obst_detect(particle_map[X][i],particle_map[Y][i])       ||
           ( (particle_map[X][i]>MAX_X) || (particle_map[X][i]<0.0) ) ||
           ( (particle_map[Y][i]>MAX_Y) || (particle_map[Y][i]<0.0) )  )
		{
	     	particle_map[W][i] = 0.0;
        }
		else
		{
	         for(j=0;j<MAX_ANGLE*2/ANGLE_STEP+1;j++)
			 {
				//s_part[j]=getSonar(particle_map[X][i],particle_map[Y][i],particle_map[Q][i]+j*ANGLE_STEP-MAX_ANGLE);
	            s_part[j]=getSonar2(particle_map[X][i],particle_map[Y][i],particle_map[Q][i]+j*ANGLE_STEP-MAX_ANGLE);
				diff = abs( sonar[j] - s_part[j] );
	            /*if( diff >=50)
	            	discrep += 50;
	            else discrep += diff;*/
				discrep += diff;
	         }
	         particle_map[W][i] = 1 - discrep/(CUTOFF_RANGE/3*(MAX_ANGLE*2/ANGLE_STEP+1));

	         discrep=0;
		 }
	}
}    

//print particle positions at each step 
void print_file(int num){
     int i;
     
     FILE *outp;
	 char name[30];
	 sprintf(name,"output/Particle_PDF%d.dat",num);
     outp = fopen(name,"w");
     
     for(i=0;i<NUM_PARTICLES;i++)
     fprintf(outp,"%.2f\t %.2f\n",particle_map[X][i],particle_map[Y][i]);
     
     fclose(outp);
}

//print robot location at each step
void print_path(int num){  
     FILE *outp;
	 char name[30];
	 sprintf(name,"output/particle%d.txt",num);
     outp = fopen(name,"w");

     fprintf(outp,"%.2f\t %.2f\n",tx,ty);
     
     fclose(outp);
}

//print initial distribution  
void print_init(void){
     int i;
     
     FILE *outp;
     outp = fopen("output/Particle_PDF.dat","w");
     
     for(i=0;i<NUM_PARTICLES;i++)
     fprintf(outp,"%.2f\t %.2f\t %.2f\n",particle_map[X][i],particle_map[Y][i],particle_map[W][i]);
     
     fclose(outp);
}

//generate obstacle file for visualization
void print_obst(void){
     int x,y,obst[50];

     FILE *outp;
     outp = fopen("output/Obstacles.dat","w");
     for(int i = 0; i < NUM_OBSTACLES; ++i)
	 {
	     for(x=obstacles[X][i]-OBS_SIZE_X/2;x<obstacles[X][i]+OBS_SIZE_X/2;x++)
	     	for(y=obstacles[Y][i]-OBS_SIZE_Y/2;y<obstacles[Y][i]+OBS_SIZE_Y/2;y++)
	     		fprintf(outp,"%d\t %d\n",x,y);
	 }          
     fclose(outp);    
}

//detect obstacles at x,y
int obst_detect(float x,float y)
{
    int obstacle=0,i;
	for(i = 0; i < NUM_OBSTACLES; ++i)
	{
		if(x <= obstacles[X][i] + OBS_SIZE_X/2 && x >= obstacles[X][i] - OBS_SIZE_X/2)
		{
			if(y <= obstacles[Y][i] + OBS_SIZE_Y/2 && y >= obstacles[Y][i] - OBS_SIZE_Y/2)
			{
				obstacle = 1;
				break;
			}
		}
	}   
    return obstacle;    
}

//compute mean and standard deviation of x,y,q
void generate_outputs(int num)
{
	int i = 0;
	float sum = 0.0;
	//Calculate the mean of the particles
	for(i = 0; i < NUM_PARTICLES; ++i)
	{
		sum+= particle_map[X][i];
	}
	x = sum/NUM_PARTICLES;

	sum = 0.0;
	for(i = 0; i < NUM_PARTICLES; ++i)
	{
		sum+= particle_map[Y][i];
	}
	y = sum/NUM_PARTICLES;

	sum = 0.0;
	for(i = 0; i < NUM_PARTICLES; ++i)
	{
		sum+= particle_map[Q][i];
	}
	q = sum/NUM_PARTICLES;

	//Calculate the sigma of the particles
	sum = 0.0;
	for(i = 0; i < NUM_PARTICLES; ++i)
	{
		sum+= particle_map[X][i]*particle_map[X][i];
	}
	sx = sqrt(sum/NUM_PARTICLES - x*x);

	sum = 0.0;
	for(i = 0; i < NUM_PARTICLES; ++i)
	{
		sum+= particle_map[Y][i]*particle_map[Y][i];
	}
	sy = sqrt(sum/NUM_PARTICLES - y*y);

	sum = 0.0;
	for(i = 0; i < NUM_PARTICLES; ++i)
	{
		sum+= particle_map[Q][i]*particle_map[Q][i];
	}
	sq = sqrt(sum/NUM_PARTICLES - q*q);

	if( (sx < CONVERGED) && (sy < CONVERGED) && (sq < CONVERGED) )
	{
		converged = 1;
	}
	else
	{
		converged = 0;
	}

	FILE *outp;
    char mean[30];
    sprintf(mean,"output/mean%d.dat",num);
    outp = fopen(mean,"w");
    fprintf(outp,"%.2f\t%.2f\t%.2f\t%.2f",x,y,sx,sy);
    fclose(outp);

	print_file(num);
}

int sign(int v)
{
	return v > 0 ? 1 : (v < 0 ? -1 : 0);
}

//unused sonar generation function
int getSonar(float x, float y, float q)
{

        int ret = MAX_RANGE;
        double a = q*PI/180; 
        for (int i = 0; i < MAX_RANGE; ++i)
        {
                if( obst_detect( (x+i*cos(a)), (y+i*sin(a)) ) )
                {
                        ret = i;
                        break;
                }
        }
        return ret;
}

//computer sonar readings at a given pose x,y,q
int getSonar2(float x, float y, float q)
{
	int ret = 20;
	int temp = 0;
    double a = q*PI/180;
	while( ret < CUTOFF_RANGE)
	{
		if( obst_detect( (x+ret*cos(a)), (y+ret*sin(a)) ) )
        {
			break;
		}
		ret += FIRST_STEP;
	}
	if(ret < CUTOFF_RANGE)
	{
		temp = ret;
		ret -= SECOND_STEP;
		while(ret > temp-FIRST_STEP)
		{
			if( !obst_detect( (x+ret*cos(a)), (y+ret*sin(a)) ) )
	        {
				break;
			}
			ret -= SECOND_STEP;
		}
		temp = ret;
		ret += 1;
		while(ret < temp+SECOND_STEP)
		{
			if( obst_detect( (x+ret*cos(a)), (y+ret*sin(a)) ) )
	        {
				break;
			}
			ret += SECOND_STEP;
		}
	}
    return ret;
}

//calculate distance between two points
float distance(float x0,float y0, float x1,float y1)
{
	float xdiff = x1 - x0;
	float ydiff = y1 - y0;
	return sqrt(xdiff*xdiff + ydiff*ydiff);
}
