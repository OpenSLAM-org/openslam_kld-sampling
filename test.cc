/*********************************************************************

  KLD-SAMPLING: Adequately Sampling from an Unknown Distribution.
  
  Copyright (C) 2006 - Patrick Beeson  (pbeeson@cs.utexas.edu)


  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
  USA

*********************************************************************/


#include <math.h>
#include <iostream>
#include <sys/time.h>
#include "kld-sampling.hh"

using namespace std;


float quantile=0.5;
float kld_error = 0.1;
float bin_size = 0.1;
int min_samples=10;
int seed=-1;
float umean=0;
float uvar=1;


// Return real between (0,1).
inline float real_random(float multi=1.0){
  return float(random())/RAND_MAX*multi;
}


// Return a point drawn from a 1D Gaussian distribution centered at
// mean with a given standard deviation.
float get_sample(float mean, float std) {
  
  // polar form of a gaussian distribution from
  // http://www.taygeta.com/random/gaussian.html

  float x1, x2, w, y1;
  static float y2; 
  static bool _ready=false;

  
  if (_ready) {
    _ready=false;
    return y2*std+mean;
  }
  
  _ready=true;
  
  do {
    x1 = 2.0 * real_random() - 1.0;
    x2 = 2.0 * real_random() - 1.0;
    w = x1 * x1 + x2 * x2;
  } while (w>1.0 || w==0.0);
  
  w = sqrt((-2.0 * log(w))/w );
  y1 = x1 * w;
  y2 = x2 * w;

  float tmp=y1*std+mean;
  
  return  tmp;
}  


// Given 1D samples, computes mean.
float get_mean(vector<float> samps) {

  float sum=0;

  for (unsigned int i=0;i<samps.size();i++)
    sum+=samps[i];

  return sum/samps.size();
}


// Given 1D samples and mean, computes variance.
float get_variance(vector<float> samps, float mean) {

  int sz=int(samps.size());

  if (sz < 2) 
    return 0;
  
  float sum=0;

  for (int i=0;i<sz;i++)
    sum+=pow(samps[i]-mean,2);
  
  return sum/(sz-1);
}


void parse_params(int argc, char *argv[]) {
  for (int i=1;i<argc;i++) {
    if (0==strcmp(argv[i],"-?")) {
      cout << "\nTo run : ./test <options>\n\n";
      cout << "options (see README for details):\n";
      cout << "-quantile Q\n";
      cout << "-error E\n";
      cout << "-bin-size B\n";
      cout << "-min-samples M\n";
      cout << "-underlying-mean U\n";
      cout << "-underlying-var V\n";
      cout << "-seed S\n";

      exit (0);
    }
    else 
      if (i+1 < argc) {
	if (0==strcmp(argv[i],"-quantile"))
	  quantile=atof(argv[++i]);
	else
	  if (0==strcmp(argv[i],"-error"))
	    kld_error=atof(argv[++i]);
	  else
	    if (0==strcmp(argv[i],"-bin-size"))
	      bin_size=atof(argv[++i]);
	    else
	      if (0==strcmp(argv[i],"-min-samples"))
		min_samples=atoi(argv[++i]);
	      else
		if (0==strcmp(argv[i],"-seed"))
		  seed=atoi(argv[++i]);
		else
		  if (0==strcmp(argv[i],"-underlying-mean"))
		    umean=atof(argv[++i]);
		  else
		    if (0==strcmp(argv[i],"-underlying-var"))
		      uvar=atof(argv[++i]);
		    else {
		      cout << "Incorrect parameter list.\n";
		      cout << "Please run with -? for runtime options.\n";
		      exit(0);
		    }
      }
      else {
	cout << "Incorrect parameter lists.\n";
	cout << "Please run with -? for runtime options.\n";
	exit(0);
      }
  }
}


int main(int argc, char *argv[]) {
  parse_params(argc,argv);

  if (quantile < 0.5 || quantile > 1) {
    cout << "quantile must be between 0.5 and 1.0\n";
    cout << "quantile is max thresholded at 0.99998\n";
    cout << "Please run with -? for runtime options.\n";
    exit (0);
  }

  quantile=fmin(0.99998,quantile);

  if (min_samples < 10) {
    cout << "min-samples needs to be at least 10.\n";
    cout << "Please run with -? for runtime options.\n";
    exit (0);
  }

  if (kld_error <= 0) {
    cout << "error must be greater than 0.\n";
    cout << "Please run with -? for runtime options.\n";
    exit (0);
  }

  if (uvar < 0) {
    cout << "underlying-var must be positive.\n";
    cout << "Please run with -? for runtime options.\n";
    exit (0);
  }

  if (bin_size <= 0) {
    cout << "bin-size must be greater than 0.\n";
    cout << "Please run with -? for runtime options.\n";
    exit (0);
  }

  if (seed == -1) {
    struct timeval tv;
    gettimeofday(&tv, NULL);    
    seed=int(tv.tv_usec);
  }
  srandom(seed);



  cout << endl<< "Source distribution: 1D Gaussian with mean=" <<umean<< " and variance=" << uvar << endl;
  cout << "KLD quantile: " << quantile << endl;
  cout << "KLD error: " << kld_error << endl;
  cout << "KLD bin size: " << bin_size << endl;
  cout << "Minimum # of samples: " << min_samples << endl;
  cout << "Random Seed: "<<seed<<endl<<endl;


  // Make into a vector of bins because the kld_sampling module
  // assumes multivariate distributions.
  vector<float> bins;
  bins.push_back(bin_size);

  kld_sampling sampler;
  sampler.init(quantile,kld_error,bins,min_samples);

  int num_samples=0;
  vector<float> samples;

  float curr_sample;
  vector<float> curr_sample2;
  curr_sample2.resize(1);

  float ustd=sqrt(uvar);

  while (num_samples < min_samples) {
    curr_sample=get_sample(umean,ustd);

    samples.push_back(curr_sample);
    num_samples++;

    //make the sample into a 1D vector because the kld_sampling module
    //assumes multivariate distributions.
    curr_sample2[0]=curr_sample;

    min_samples=sampler.update(curr_sample2);
  }

  float mean=get_mean(samples);
  float variance=get_variance(samples,mean);

  cout << "Final number of samples: "<<num_samples<<endl;
  cout << "Final mean: "<<mean<<endl;
  cout << "Final variance: "<<variance<<endl<<endl;

  return 1;
}
