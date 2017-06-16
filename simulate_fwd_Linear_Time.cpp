# include <vector>
# include <iostream>
//# include <random>
# include <algorithm>
# include <cmath>
# include <math.h>
# include <stdlib.h>
# include <fstream>

using namespace std; 

class particle
{
	public:
		double x;
		double y;
		double interference ;
		particle(double lx,double ly)
		{
			x = lx;
			y = ly;
			interference = 0;

		}
};

class system_class // The class containing the system of particles
{
public:
	

	double S ;
	double L ;
	std::vector<particle *> particles;
	std::vector<double> rates;
	std::vector<double> death_probability ;
	std::vector<double> cumulative_death_probability;
	double total_rate;
	double noise ;
	double lambda ;
	int death_count;
	int birth_count;
	
	


	double get_uniform_random_number() // Returns a random number to be used by the program.
	{
		return (double)rand()/RAND_MAX;
		
	}

	
	system_class(int initial_points) // Initialize the System 
	{

		// System Parameters - PLAY WITH
		noise = 0.1;
		lambda = 0.47;
		S = 5;
		L = 1; 
		death_count = 0;
		birth_count = 0;


		// Initializing the random seed 
		time_t beginning_time;
		time(&beginning_time); 

		srand(time(NULL));



		particles.resize(0);
		rates.resize(0);
		death_probability.resize(0);
		cumulative_death_probability.resize(0);

		
		for(int loopvar = 0; loopvar < initial_points; loopvar++)
		{
			double lx = get_uniform_random_number();
			double ly = get_uniform_random_number();
			particle *local_particle = new particle(S*lx,S*ly);
			particles.push_back(local_particle);

		}
		for(int lvar = 0; lvar < particles.size();lvar++)
		{
			double linterference = 0;
			for(int l2var = 0; l2var < particles.size(); l2var++)
			{
				if(l2var == lvar)
					continue;
				double xd1 = std::min(std::abs(particles[lvar]->x - particles[l2var]->x) , S - std::abs(particles[lvar]->x - particles[l2var]->x) );
				double yd1 = std::min(std::abs(particles[lvar]->y - particles[l2var]->y) , S - std::abs(particles[lvar]->y - particles[l2var]->y) );
				linterference += path_loss(pow(  pow(xd1,2)+pow(yd1,2)  , 0.5));

			}
			particles[lvar]->interference = linterference;
		}

	}

	

	double path_loss(double r)
	{
		return pow( r+1, -2.25);  // For S=5, the stability boundary is at 0.84
	}

	void compute_probability() // Compute the death probability vector at a time step
	{
		
		int N = particles.size();
		rates.resize(0);
		death_probability.resize(0);
		cumulative_death_probability.resize(0);
		total_rate = 0.0;

		//THE INTERFERNCE IN A CONFIGURATION IS COMPUTED ALREADY

		
		
		double run_sum = 0;

		for(int loopvar1 = 0;loopvar1 < N; loopvar1++)
			total_rate += log10(1 + (1/(noise + particles[loopvar1]->interference)))/log10(2) ;


		for(int loopvar = 0; loopvar < N; loopvar++)
		{
			//run_sum += (rates[loopvar]/total_rate);
			run_sum += (log10(1 + (1/(noise + particles[loopvar]->interference)))/log10(2))/total_rate ;
			death_probability.push_back((log10(1 + (1/(noise + particles[loopvar]->interference)))/log10(2))/total_rate );
			cumulative_death_probability.push_back(run_sum);

			
		}
	
	}


	void generate_next_time_step()
	{
		// Figure out if a birth or death event happens and generate the next particle vector.
		double total_birth_rate = lambda*S*S;
		
		double total_death_rate = total_rate ;
		
		double total_event_rate = total_birth_rate + total_death_rate ;
		double U = get_uniform_random_number();
		if(U <= total_birth_rate/total_event_rate)
		{
			// Birth Happens - EASY
			

			double lx = get_uniform_random_number();
			double ly = get_uniform_random_number();
			particle *temp_particle = new particle(S*lx,S*ly);
			particles.push_back(temp_particle);
			birth_count++;
			// Re compute the interference
			double run_sum_l = 0;
			for(int lvar = 0; lvar < particles.size()-1; lvar++)
			{
				double xdist = std::min(std::abs(particles[lvar]->x - S*lx) , S - std::abs(particles[lvar]->x - S*lx) );
				double ydist = std::min(std::abs(particles[lvar]->y - S*ly) , S - std::abs(particles[lvar]->y - S*ly) );
				particles[lvar]->interference += path_loss(pow(  pow(xdist,2) + pow(ydist,2) , 0.5));
				run_sum_l += path_loss(pow(  pow(xdist,2) + pow(ydist,2) , 0.5));
			}
			particles[particles.size()-1]->interference = run_sum_l;

		}
		else
		{

			// Get Index of the Killed point and Remove It
			//cout<<" Death Happens with U value "<<U<<endl;
			double UU = get_uniform_random_number();
			int Nppl = particles.size();
			int loopvar;
			for(loopvar = 0; loopvar < Nppl; loopvar++)
			{
				//cout<<cumulative_death_probability[loopvar]<<" "<<loopvar+1<<" UU "<<UU<<" Index "<<loopvar<<endl;
				if (UU <= cumulative_death_probability[loopvar]) // loopvar is the index killed
				{
					// cout<<" UU value "<<UU<<endl;
					// cout<<" The after value "<<cumulative_death_probability[loopvar]<<endl;
					// cout<<"Dead Location Before "<<particles[loopvar]->x<<" , "<<particles[loopvar]->y<<endl;
					double xdist_d = particles[loopvar]->x;
					double ydist_d = particles[loopvar]->y;

					particles[loopvar] = particles.back();
					// Recompute the interference

					particles.pop_back();
					for (int lvar = 0; lvar < particles.size(); lvar++)
					{
						double xd1 = std::min(std::abs(particles[lvar]->x - xdist_d) , S - std::abs(particles[lvar]->x - xdist_d) );
						double yd1 = std::min(std::abs(particles[lvar]->y - ydist_d) , S - std::abs(particles[lvar]->y - ydist_d) );
						particles[lvar]->interference -= path_loss(pow(  pow(xd1,2)+pow(yd1,2)  , 0.5));
					}
					death_count--;
					break;
				}
				

			}
			
			
		}

	
	}
	void run_system(int max_time)
	{
		// Run the system forward from time 0 till time max_time
		for(int timestep = 1; timestep <= max_time; timestep++)
		{
			compute_probability();
			generate_next_time_step();
			if(timestep%500 == 0)
			   cout<<"Finished Time Step  "<<timestep<<"  with particle size "<<particles.size()<<endl;
		}
		cout<<"The size of the population at end is   "<<particles.size()<<endl;
		cout<<" The Death Count is "<<death_count<<endl;
		cout<<" The birth Cout is "<<birth_count<<endl;
		collect_statistics();


	}

	void collect_statistics()
	{
		
		ofstream MyFile;
		MyFile.open("Particle_Locations_47.csv",ios::out | ios::ate | ios::app | ios:: binary);
		for(int i = 0;i < particles.size();i++)
		{
			MyFile<<particles[i]->x<<", "<<particles[i]->y<<endl;
		}
		MyFile.close();



	}
	



};

int main(void)
{
	
	system_class *system = new system_class(0) ;
	system->run_system(45500000);
	

	

	//cout << system->get_uniform_random_number() << endl;
	//cout << system->get_uniform_random_number() << endl;


}
