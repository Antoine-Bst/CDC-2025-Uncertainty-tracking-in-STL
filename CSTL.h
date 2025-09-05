/* ============================================================================
 * D Y N I B E X - STL formula verification on tube
 * ============================================================================
 * Copyright   : ENSTA
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Antoine Besset, Joris Tillet and Julien Alexandre dit Sandretto
 * Created     : Jul 22, 2025
 * Modified    : Jul 22, 2025
 * Sponsored   : This research benefited from the support of the "STARTS Projects - CIEDS - Institut Polytechnique"
 * ---------------------------------------------------------------------------- */
#ifndef CSTL_H
#define CSTL_H

#include <iostream>
#include <vector>
#include "ibex.h"
#include <unordered_set>

using namespace std;
using namespace ibex;

typedef pair<pair<int, vector<int>>, pair<double, double>> TimeInterval;

// Function prototypes
vector<TimeInterval> neg_stl_titv(const vector<TimeInterval>& sp);
vector<TimeInterval> and_stl_titv(const vector<TimeInterval>& sp1, const vector<TimeInterval>& sp2);
vector<TimeInterval> or_stl_titv(const vector<TimeInterval>& sp1, const vector<TimeInterval>& sp2);
void print_TimeIntervals(const vector<TimeInterval>& TimeIntervals);
vector<TimeInterval> until_stl_minkowski(const vector<TimeInterval>& list1, const vector<TimeInterval>& list2, pair<double, double> time_itv);
vector<TimeInterval> Finally(const vector<TimeInterval>& list1, pair<double, double> time_itv);
vector<TimeInterval> Globally(const vector<TimeInterval>& list1, pair<double, double> time_itv);
vector<int> minusonecleaner(const vector<int>& input);
vector<TimeInterval> merge_TimeIntervals_Mink(const vector<TimeInterval>& TimeIntervals);
vector<TimeInterval> merge_TimeIntervals(const vector<TimeInterval>& TimeIntervals);
void print_TimeIntervals_to_file(const vector<TimeInterval>& TimeIntervals, const string& file_path);

//tools
std::vector<std::pair<IntervalVector, Interval>> Sim_to_jn_tube(ibex::simulation& sim);
std::pair<std::vector<vector<TimeInterval>>, std::vector<Interval>> predicate_satisfaction_jn(const std::vector<std::pair<IntervalVector, Interval>>& jn_box, const vector<IntervalVector>& predicate_list);
std::vector<std::pair<IntervalVector, Interval>> Append_tube(const std::vector<std::pair<IntervalVector, Interval>>& old_tube, const std::vector<std::pair<IntervalVector, Interval>>& extension_tube);
double last_time_tube(const std::vector<std::pair<IntervalVector, Interval>>&  result);
double first_time_tube(const std::vector<std::pair<IntervalVector, Interval>>&  result);
std::pair<std::vector<vector<TimeInterval>>, std::vector<Interval>> predicate_satisfaction(ibex::simulation& sim, const vector<IntervalVector>& predicate_list);
vector<double> Bisection_time(const vector<Interval>& Un_origin, const double& min_step_time, const double& t_min);
vector<Interval> Uncertainty_origin(vector<TimeInterval> Phi1, const vector<Interval>& P_Satisfaction_signals);
void process_and_save_data(ibex::simulation& sim, 
                           const vector<IntervalVector>& predicate_list, 
                           const vector<int>& selected_indices,
                           const string& predicate_file,
                           const string& jn_box_file);

int satisfies_at_time(const double& time, const vector<TimeInterval>& phi);
std::vector<std::pair<ibex::IntervalVector, ibex::Interval>> Sim_to_jn_tube_test(const ibex::simulation& sim);
#endif // CSTL_H

