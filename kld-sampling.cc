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


#include "kld-sampling.hh"
#include <fstream>
#include <iostream>
#include <math.h>

vector<float> kld_sampling::ztable;

// Constructs Z-table (from ztable.data) to lookup statistics.
kld_sampling::kld_sampling() {
  if (kld_sampling::ztable.empty())
    build_table();
}


/**
   Initialize a round of KLD sampling.  Takes in kld-parameters:
   quantile, kld-error, bin size, minimum number of samples
**/
void kld_sampling::init(float quantile, float err, const vector<float>& bsz, int sample_min) {
  support_samples=0;
  num_samples=0;
  if (sample_min < absolute_min)
    kld_samples=absolute_min;
  else kld_samples=sample_min;

  bins.clear();

  confidence=quantile-0.5; // ztable is from right side of mean
  confidence=fmin(0.49998,fmax(0,confidence));

  max_error=err;
  bin_size=bsz;

  zvalue=4.1;
  for (unsigned int i=0; i<kld_sampling::ztable.size();i++)
    if (kld_sampling::ztable[i] >= confidence) {
      zvalue=i/100.0;
      break;
    }
}


/**
   Update kld-sampler with the last sample drawn.  Returns a guess at
   the number of samples needed before the distribution (which is
   unknown) is adequately sampled.
**/
int kld_sampling::update(const vector<float>& sample) {
  if (bin_size.empty()) {
    cerr << "kld-sampling.cc: Must run init() before update()\n";
    exit (-1);
  }

  if (sample.size() != bin_size.size()){
    cerr << "kld-sampling.cc: Sample size not the same number of dimensions as the bins\n";
    exit(-1);
  }
    
  curr_sample=sample;
  num_samples++;

  if (in_empty_bin()) {
    support_samples++;
    if (support_samples >=2) {
      int k=support_samples-1;
      k=(int)ceil(k/(2*max_error)*pow(1-2/(9.0*k)+sqrt(2/(9.0*k))*zvalue,3));
      if (k > kld_samples)
	kld_samples=k;
    }
  }
  return kld_samples;
}


/**
   Builds a z-table which is necessary for the statiscal kld-sampling.
**/
void kld_sampling::build_table() {
  float tmp;
  ifstream ifile("ztable.data");

  if (ifile.is_open()) {
    while (!ifile.eof()) {
      ifile >> tmp;
      kld_sampling::ztable.push_back(tmp);
    }
  }
  else {
    cerr << "kld-sampling.cc: ztable.data does not exist. Error!\n";
    exit(-1);
  }
  
}


/**
   Determines whether a sample falls into a bin that has already been
   sampled.
 **/
bool kld_sampling::in_empty_bin() {

  vector<float> curr_bin;
  
  for (unsigned int i=0; i<curr_sample.size(); i++)
    curr_bin.push_back(floor(curr_sample[i]/bin_size[i]));

  for (unsigned int i=0; i < bins.size(); i++)
    if (curr_bin==bins[i]) {
      return false;
    }
  bins.push_back(curr_bin);
  return true;
}
