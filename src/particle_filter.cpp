/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	num_particles = 100;
	std::default_random_engine gen;
	
	std::normal_distribution<double> N_x(x, std[0]);
	std::normal_distribution<double> N_y(y, std[1]);
	std::normal_distribution<double> N_theta(theta, std[2]);
	
	for (int i = 0; i < num_particles; i++) {
		Particle particle;
		particle.id = i;
		particle.x = N_x(gen);
		particle.y = N_y(gen);
		particle.theta = N_theta(gen);
		particle.weight = 1;

		particles.push_back(particle);
		weights.push_back(particle.weight);
	}
	
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	default_random_engine gen;
	
	for (int i = 0; i < num_particles; i++) 
	{
		double new_y;
		double new_x;
		double new_theta;

		if (yaw_rate == 0) 
		{
			new_x = particles[i].x + velocity*delta_t*cos(particles[i].theta);
			new_y = particles[i].y + velocity*delta_t*sin(particles[i].theta);
			new_theta = particles[i].theta;
		} 
		else 
		{
			new_x = particles[i].x + (velocity/yaw_rate)*(sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
			new_y = particles[i].y + (velocity/yaw_rate)*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
			new_theta = particles[i].theta + yaw_rate*delta_t;
		}

		normal_distribution<double> N_x(new_x, std_pos[0]);
		normal_distribution<double> N_y(new_y, std_pos[1]);
		normal_distribution<double> N_theta(new_theta, std_pos[2]);

		particles[i].x = N_x(gen);
		particles[i].y = N_y(gen);
		particles[i].theta = N_theta(gen);
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	for(int i = 0; i < observations.size(); i++)
	{
		// Set minimum distance to biggest value possible
		double minimum_distance = numeric_limits<double>::max();

		for(int j = 0; j < predicted.size(); j++)
		{
			//double diffX = observations[i].x - predicted[j].x;
			//double diffY = observations[i].y - predicted[j].y;
			//double current_distance = sqrt(pow(diffX,2)+pow(diffY,2));
			
			double current_distance = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
			if(current_distance < minimum_distance)
			{
				minimum_distance = current_distance;
				observations[i].id = predicted[j].id;
			}
		}
	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) 
{
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	

	/*1
	Need to:
		1. Find landmarks
		2. Transform observation coordinates
		3. associate observations to landmark
		4. Calculate weights
	*/
	for(int p = 0; p < num_particles; p++)
	{
		vector<double> sense_x;
		vector<double> sense_y;
		
		// Find landmarks
		vector<LandmarkObs> landmarks_found;
		for(int i = 0; i < map_landmarks.landmark_list.size(); i++)
		{
			float x_landmarks = map_landmarks.landmark_list[i].x_f;
			float y_landmarks = map_landmarks.landmark_list[i].y_f;
			int id_landmarks = map_landmarks.landmark_list[i].id_i;

			// Use pythagoras theorem for checking sensor range
			float diff_x2 = pow(particles[p].x - x_landmarks, 2);
			float diff_y2 = pow(particles[p].y - y_landmarks, 2);

			//if((diff_x2 + diff_y2) <= pow(sensor_range, 2))
			if(fabs(x_landmarks - particles[p].x)<= sensor_range && fabs(y_landmarks - particles[p].y)<= sensor_range)
				landmarks_found.push_back(LandmarkObs{id_landmarks, x_landmarks, y_landmarks});

		}


		// Transform obseration coordinates
		vector<LandmarkObs> trans_observations;
		LandmarkObs obs;
		for(int i = 0; i < observations.size(); i++)
		{
			LandmarkObs trans_obs;
			
			// Perform the space transformation from vehicle tomap
			trans_obs.x = particles[p].x + (observations[i].x*cos(particles[p].theta) - observations[i].y*sin(particles[p].theta));
			trans_obs.y = particles[p].y + (observations[i].x*sin(particles[p].theta) + observations[i].y*cos(particles[p].theta));
			trans_observations.push_back(trans_obs);
		}
		
		// Associate observations to landmark
		dataAssociation(landmarks_found, trans_observations);
		

		// Calculate weights
		particles[p].weight = 1;
		for(int i = 0; i < trans_observations.size(); i++)
		{
			double meas_x = trans_observations[i].x;	// Transformed observation x landmark
			double meas_y = trans_observations[i].y;	// Transformed observation y landmark
			double std_x = std_landmark[0];
			double std_y = std_landmark[1];

			for(int j = 0; j < landmarks_found.size(); j++)
			{
				double mu_x = map_landmarks.landmark_list[j].x_f;	// Landmark x that was located
				double mu_y = map_landmarks.landmark_list[j].y_f;	// Landmark y that was located

				if (landmarks_found[j].id == trans_observations[i].id)
				{
					long double dx2 = pow(meas_x-mu_x,2);
					long double dy2 = pow(meas_y-mu_y,2);

					long double std_x2 = pow(std_x,2);
					long double std_y2 = pow(std_y,2);

					long double exp_equation = dx2/(2*std_x2) + dy2/(2*std_y2);

					long double multipler = 1/(2*M_PI*std_x*std_y)*exp(-exp_equation);
					
					if(multipler > 0)
						particles[p].weight *= multipler;
					
					//associations.push_back(association+1);
					//sense_x.push_back(trans_observations[i].x);
					//sense_y.push_back(trans_observations[i].y);
				}
				weights[p] = particles[p].weight;
			}
			
		}

		// Normalize weights
		long double weights_sum = 0;
		for(int i = 0; i < weights.size(); i++)
			weights_sum += weights[p];
		
		for(int i = 0; i < weights.size(); i++)
			weights[p] /= weights_sum;
		
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	default_random_engine gen;
	discrete_distribution<int> distribution(weights.begin(), weights.end());
	
	vector<Particle> resample_particles;
	
	for(int i = 0; i < num_particles; i++)
		resample_particles.push_back(particles[distribution(gen)]);
	
	particles = resample_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

	particle.associations.clear();
    particle.sense_x.clear();
    particle.sense_y.clear();
	
    particle.associations = associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
	
	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
