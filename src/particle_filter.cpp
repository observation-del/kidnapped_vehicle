/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles         = 20;  // TODO: Set the number of particles
  particles.resize(num_particles);

  std::default_random_engine gen;
  double std_x, std_y, std_theta;  // Standard deviations for x, y, and theta

  // Set standard deviations for x, y, and theta
  std_x                 = std[0];
  std_y                 = std[1];
  std_theta             = std[2];
  
  // Create normal distributions for y and theta
  std::normal_distribution<double> dist_x(0, std_x);
  std::normal_distribution<double> dist_y(0, std_y);
  std::normal_distribution<double> dist_theta(0, std_theta);
  unsigned int i = 0;
  for(auto &e : particles){
    e.id              = i;
    e.x               = x + dist_x(gen);
    e.y               = y + dist_y(gen);
    e.theta           = theta + dist_theta(gen);
    e.weight          = 1.0;
    i++;
  }

  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

  std::default_random_engine gen;
  double std_x, std_y, std_theta;  // Standard deviations for x, y, and theta
  double dx, dy, dtheta;
  // Set standard deviations for x, y, and theta
  std_x                 = std_pos[0];
  std_y                 = std_pos[1];
  std_theta             = std_pos[2];
  

  for(auto &e : particles){
    if(yaw_rate != 0){
      dx                  = (velocity/yaw_rate) * (sin(e.theta+yaw_rate*delta_t) - sin(e.theta));
      dy                  = (velocity/yaw_rate) * (cos(e.theta)-cos(e.theta+yaw_rate*delta_t));
      dtheta              = yaw_rate * delta_t;
    }
    
    else{
      dx                  = velocity * delta_t * cos(e.theta);
      dy                  = velocity * delta_t * sin(e.theta);
      dtheta              = 0.0;
    }
    
    // Create normal distributions for y and theta
    std::normal_distribution<double> dist_x(dx, std_x);
    std::normal_distribution<double> dist_y(dy, std_y);
    std::normal_distribution<double> dist_theta(dtheta, std_theta);

    e.x                 += dist_x(gen);
    e.y                 += dist_y(gen);
    e.theta             += dist_theta(gen);

  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

  double dist;
  double tmp_dist;

  for(auto &observation : observations){
    dist = INFINITY;
    for(auto &e : predicted){
      tmp_dist = std::pow((observation.x - e.x),2.0) + std::pow((observation.y - e.y), 2.0);
      if(tmp_dist < dist){
        observation.id      = e.id;
        dist                = tmp_dist;
      }
    }
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  std::vector<LandmarkObs> predictions;
  vector<LandmarkObs> observations_T;
  LandmarkObs observation_T;
  double sum_weight = 0.0;
  weights.resize(0);
  for(auto &particle : particles){
    particle.weight = 1.0;
    // reset predictions
    predictions.resize(0);
    // Each map landmark for loop
    for(auto &landmark : map_landmarks.landmark_list){
      //Get id and x,y coordinates
      float lm_x = landmark.x_f;
      float lm_y = landmark.y_f;
      int lm_id = landmark.id_i;
      if(fabs(lm_x - particle.x) <= sensor_range && fabs(lm_y - particle.y) <= sensor_range){
        predictions.emplace_back(LandmarkObs{ lm_id, lm_x, lm_y});
      }
    }

    observations_T.resize(0);
    for(auto &observation : observations){
    // rotation
    observation_T.x = observation.x * cos(particle.theta) - observation.y * sin(particle.theta);
    observation_T.y = observation.x * sin(particle.theta) + observation.y * cos(particle.theta);

    // translation
    observation_T.x += particle.x;
    observation_T.y += particle.y;

    observations_T.emplace_back(observation_T);
    }

    // assign id
    if(predictions.size() > 0){
      dataAssociation(predictions, observations_T);
      for(auto &observation : observations_T){
        std::vector<LandmarkObs>::iterator itr;
        std::size_t index;
        itr = std::find_if(predictions.begin(), predictions.end(), [&](LandmarkObs check) { return check.id == observation.id; });
        index = std::distance(predictions.begin(), itr);
        double weight = multiv_prob(std_landmark[0], std_landmark[1], observation.x, 
          observation.y, predictions[index].x, predictions[index].y);
        particle.weight *= weight;
      }
      sum_weight += particle.weight;
    }
  } 

  // normalize
  for(auto &particle : particles){
    particle.weight /= sum_weight;
    weights.emplace_back(particle.weight);
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  std::default_random_engine gen;
  std::discrete_distribution<size_t> index(weights.begin(), weights.end());

  // Set of new particles
  std::vector<Particle> new_particles;

  for(auto &particle : particles){
    new_particles.emplace_back(particles[index(gen)]);
  }
  particles.resize(0);
  std::copy(new_particles.begin(), new_particles.end(), std::back_inserter(particles));
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}