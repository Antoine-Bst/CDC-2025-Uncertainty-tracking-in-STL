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

#include <iostream>
#include <vector>
#include "ibex.h"
#include "CSTL.h"
#include <unordered_set>

using namespace std;
using namespace ibex;

typedef pair<pair<int, vector<int>>, pair<double, double>> TimeInterval;
//typedef vector<vector<int>> Undef_index;

//Undef_index Tableofunkown;
//on marque les index des incertitudes au début et à chaque fois on renvoie une liste de branche
//en gros pour chaque incertitudes indexer je dois connaitre au niveau du predicat d'où vient l'incertitude. ça me permet de bisecter par la suite juste la boite ou ça merde
//quand y'a des ou ou des until je dois être capable de virer des index

int and_unitary(int val1, int val2) {
    return (val1 > 0 && val2 > 0) ? max(val1, val2) : 0;
}

int or_unitary(int val1, int val2) {
    if (val1 > 0 && val2 > 0) return min(val1, val2);
    if (val1 == 0 && val2 > 0) return val2;
    if (val2 == 0 && val1 > 0) return val1;
    return 0;
}

int neg_unitary(int val){
    if (val==1)
    {
        return 0;
    }
    else if (val==0)
    {
        return 1;
    }
    else
    {
        return 2;
    } 
}

std::vector<int> concatVectors(const std::vector<int>& vec1, const std::vector<int>& vec2) {
    std::vector<int> result;
    result.reserve(vec1.size() + vec2.size()); // Optimisation pour éviter des reallocations

    result.insert(result.end(), vec1.begin(), vec1.end()); // Ajouter vec1
    result.insert(result.end(), vec2.begin(), vec2.end()); // Ajouter vec2

    return result;
}

std::vector<int> minusonecleaner(const std::vector<int>& input) {
    std::unordered_set<int> seen;  // To track unique values
    std::vector<int> temp;

    for (int num : input) {
        if (num != -1 && seen.find(num) == seen.end()) { 
            seen.insert(num); // Mark as seen
            temp.push_back(num);
        }
    }

    return temp;
}

vector<TimeInterval> completator(const vector<TimeInterval>& sp, double fin_max) {
    if (sp.empty()) return {{{2,{-2}}, make_pair(0, fin_max)}};
    vector<TimeInterval> resultat = sp;
    if (sp.back().second.second<fin_max)
    {
        resultat.push_back({{2,{-2}}, make_pair(sp.back().second.second, fin_max)}); //unmarked uncertainty (no knowledge)
    }
    return resultat;
}

vector<TimeInterval> merge_TimeIntervals(const vector<TimeInterval>& TimeIntervals) {
    if (TimeIntervals.empty()) return {};

    vector<TimeInterval> merged;
    merged.push_back(TimeIntervals[0]);

    for (size_t i = 1; i < TimeIntervals.size(); ++i) {
        if (TimeIntervals[i].first.first == merged.back().first.first) {
            merged.back().second.second = TimeIntervals[i].second.second;
            if (TimeIntervals[i].first.first == 2)
            {
                merged.back().first.second = concatVectors(TimeIntervals[i].first.second, merged.back().first.second);
            }
            else
            {
                merged.back().first.second = {-1};
            }
            
            
        } else {
            merged.push_back(TimeIntervals[i]);
        }
    }
    //std::cout<<"jemerge"<<std::endl;
    return merged;
}


vector<TimeInterval> merge_TimeIntervals_Mink(const vector<TimeInterval>& TimeIntervals) {
    if (TimeIntervals.empty()) return {};

    vector<TimeInterval> merged;
    merged.push_back(TimeIntervals[0]);

    for (size_t i = 1; i < TimeIntervals.size(); ++i) {
        if (TimeIntervals[i].first.first>0 && merged.back().first.first>0) {
            merged.back().second.second = TimeIntervals[i].second.second;
            merged.back().first.first = 2;
            merged.back().first.second = concatVectors(TimeIntervals[i].first.second, merged.back().first.second);
        } else {
            merged.push_back(TimeIntervals[i]);
        }
    }
    
    return merged;
}

vector<TimeInterval> neg_stl_titv(const vector<TimeInterval>& sp) {
    vector<TimeInterval> negated_TimeIntervals;
    for (size_t i = 0; i < sp.size(); i++)
    {
       negated_TimeIntervals.push_back({{neg_unitary(sp[i].first.first),{sp[i].first.second}}, sp[i].second});
    }
    return negated_TimeIntervals;
}

vector<TimeInterval> and_stl_titv(const vector<TimeInterval>& sp1, const vector<TimeInterval>& sp2) {
    double fin_max = max(sp1.empty() ? 0 : sp1.back().second.second, sp2.empty() ? 0 : sp2.back().second.second);
    vector<TimeInterval> sp3 = completator(sp1, fin_max);
    vector<TimeInterval> sp4 = completator(sp2, fin_max);
    vector<TimeInterval> resultat;
    size_t i = 0, j = 0;

    while (i < sp3.size() && j < sp4.size()) {
        int v1 = sp3[i].first.first;
        int v2 = sp4[j].first.first;
        double debut_inter = max(sp3[i].second.first, sp4[j].second.first);
        double fin_inter = min(sp3[i].second.second, sp4[j].second.second);
        
        if (debut_inter < fin_inter) {
            if (and_unitary(v1, v2)==2)
            {
                resultat.push_back({{and_unitary(v1, v2),concatVectors(sp3[i].first.second, sp4[j].first.second)}, {debut_inter, fin_inter}});
            }
            else
            {
                resultat.push_back({{and_unitary(v1, v2),{-1}}, make_pair(debut_inter, fin_inter)});
            }
        }
        
        if (sp3[i].second.second <= sp4[j].second.second) ++i;
        else ++j;
    }
    return merge_TimeIntervals(resultat);
}


vector<TimeInterval> or_stl_titv(const vector<TimeInterval>& sp1, const vector<TimeInterval>& sp2) {
    double fin_max = max(sp1.empty() ? 0 : sp1.back().second.second, sp2.empty() ? 0 : sp2.back().second.second);
    vector<TimeInterval> sp3 = completator(sp1, fin_max);
    vector<TimeInterval> sp4 = completator(sp2, fin_max);
    vector<TimeInterval> resultat;
    size_t i = 0, j = 0;

    while (i < sp3.size() && j < sp4.size()) {
        int v1 = sp3[i].first.first;
        int v2 = sp4[j].first.first;
        double debut_inter = max(sp3[i].second.first, sp4[j].second.first);
        double fin_inter = min(sp3[i].second.second, sp4[j].second.second);
        
        if (debut_inter < fin_inter) {
            if (or_unitary(v1, v2)==2)
            {
                resultat.push_back({{or_unitary(v1, v2),concatVectors(sp3[i].first.second, sp4[j].first.second)}, {debut_inter, fin_inter}});
            }
            else
            {
                resultat.push_back({{or_unitary(v1, v2),{-1}}, make_pair(debut_inter, fin_inter)});
            }
            
        }
        
        if (sp3[i].second.second <= sp4[j].second.second) ++i;
        else ++j;
    }
    return merge_TimeIntervals(resultat);
}

vector<TimeInterval> offset_time(const vector<TimeInterval>& intervals, double offset) {
    vector<TimeInterval> offset_intervals;
    for (const auto& interval : intervals) {
        offset_intervals.push_back(make_pair(interval.first, make_pair(interval.second.first + offset, interval.second.second + offset)));
    }
    return offset_intervals;
}

vector<TimeInterval> merge_stl_titv(const vector<TimeInterval>& sp1, const vector<TimeInterval>& sp2) {
    double fin_max = max(sp1.empty() ? 0 : sp1.back().second.second, sp2.empty() ? 0 : sp2.back().second.second);
    vector<TimeInterval> sp3 = completator(sp1, fin_max);
    vector<TimeInterval> sp4 = completator(sp2, fin_max);

    vector<TimeInterval> merged_result;
    size_t i = 0, j = 0;

    while (i < sp3.size() && j < sp4.size()) {
        int v1 = sp3[i].first.first;
        int v2 = sp4[j].first.first;
        double d1 = sp3[i].second.first, f1 = sp3[i].second.second;
        double d2 = sp4[j].second.first, f2 = sp4[j].second.second;

        double start_inter = max(d1, d2);
        double end_inter = min(f1, f2);

        if (start_inter < end_inter) {
            int merged_value;
            if (v1 == 2) {
                merged_value = v2;
            } else if (v2 == 2) {
                merged_value = v1;
            } else {
                merged_value = v1;
            }

            if (merged_value == 2)
            {
            merged_result.push_back({{merged_value, concatVectors(sp3[i].first.second, sp4[j].first.second)}, make_pair(start_inter, end_inter)});
            }
            else
            {
            merged_result.push_back({{merged_value, {-1} }, make_pair(start_inter, end_inter)});
            }
                        

        }

        if (f1 <= f2) {
            i++;
        } else {
            j++;
        }
    }
    return merge_TimeIntervals(merged_result);
}

void print_TimeIntervals(const vector<TimeInterval>& TimeIntervals) {
    for (const auto& TimeInterval : TimeIntervals) {
        cout << "[" << TimeInterval.first.first << ", (" << TimeInterval.second.first << ", " << TimeInterval.second.second << ")]\n";
    }
}

void print_TimeIntervals_to_file(const vector<TimeInterval>& TimeIntervals, const string& file_path) {
    ofstream outFile(file_path); // Open file for writing

    if (!outFile) { // Check if the file is opened successfully
        cerr << "Error: Unable to open file " << file_path << endl;
        return;
    }

    for (const auto& TimeInterval : TimeIntervals) {
        outFile << "[" << TimeInterval.first.first << ", (" << TimeInterval.second.first << ", " << TimeInterval.second.second << ")]\n";
    }

    outFile.close(); // Close the file
    std::cout << "Output written to " << file_path << endl;
}


TimeInterval computeIntersection(const TimeInterval& elem1, const TimeInterval& elem2, const int& value) {
    TimeInterval result;
            double start = std::max(elem1.second.first, elem2.second.first);
            double end = std::min(elem1.second.second, elem2.second.second);
            
            if (0 <= start && start < end) { // Valid intersection
                if (elem1.first.first == 2 || elem2.first.first ==2)
                {
                   result = {{value, concatVectors(elem1.first.second, elem2.first.second)}, {start, end}};

                }
                else
                {
                    result = {{value, {-1}}, {start, end}};
                }
                //std::cout<< start <<" , " << end << std::endl;
            }
            else if (0<end && start < end) // on travaille que sur des temps positifs anyway
            {
                if (elem1.first.first == 2 || elem2.first.first ==2)
                {
                   result = {{value, concatVectors(elem1.first.second, elem2.first.second)}, {0, end}};
                }
                else
                {
                    result = {{value, {-1}}, {0, end}};
                }
            }
            else
            {
                result = {{-1, {-1}}, {-1, 0}};
            }
        return result;
}


vector<TimeInterval> until_stl_minkowski(const vector<TimeInterval>& list1, const vector<TimeInterval>& list2, pair<double, double> time_itv) {
const double fin_max = max(list1.back().second.second, list2.back().second.second);
    vector<TimeInterval> sp1 = completator(merge_TimeIntervals(list1), fin_max);
    vector<TimeInterval> sp2 = completator(merge_TimeIntervals(list2), fin_max);

double maxtime = max(sp2.back().second.second - time_itv.second, sp1.back().second.second - time_itv.second);
   std::vector<TimeInterval> result = {{{0,{-1}}, {0,maxtime}}};
   std::vector<TimeInterval> result_final = {{{0,{-1}}, {0,maxtime}}};
   TimeInterval intersection;
   std::vector<TimeInterval> temp;
   std::vector<TimeInterval> temp_final;
   std::vector<TimeInterval> merged_sp1= merge_TimeIntervals_Mink(sp1);
    
    for (size_t i = 0; i < merged_sp1.size(); i++)
    {
        for (size_t j = 0; j < sp2.size(); j++)
        {
            if (merged_sp1[i].first.first==1 && sp2[j].first.first==1){
                TimeInterval intersection = computeIntersection(merged_sp1[i], sp2[j], 1); //changer la value
                intersection.second.first = intersection.second.first-time_itv.second;
                intersection.second.second = intersection.second.second-time_itv.first;
                intersection = computeIntersection(merged_sp1[i], intersection, 1);
                
                if (intersection.first.first != -1)
                {
                temp = {{{0,{-1}},{0,intersection.second.first}}, intersection, {{0,{-1}}, {intersection.second.second,maxtime}}};
                //print_TimeIntervals(temp);
                //std::cout<<"----ici----"<<std::endl;
                result = or_stl_titv(result, temp);
                }            
            }
            else if (merged_sp1[i].first.first>0 && sp2[j].first.first>0)
            {   
                TimeInterval intersection = computeIntersection(merged_sp1[i], sp2[j], 2); //changer la value
                intersection.second.first = intersection.second.first-time_itv.second;
                intersection.second.second = intersection.second.second-time_itv.first;
                intersection = computeIntersection(merged_sp1[i], intersection, 2);
                
                if (intersection.first.first != -1)
                {   
                    //vector<int> temp_index = merged_sp1[i].first.second; //même si sp2 est à -1 on s'en fout on le vire à la fin
                    //temp_index.insert(temp_index.end(),sp2[j].first.second.begin(),sp2[j].first.second.end());
                    //intersection.first.second= temp_index;
                    //std::cout<<"----ici----"<<std::endl;
                temp = {{{0,{-1}},{0,intersection.second.first}}, intersection, {{0,{-1}}, {intersection.second.second,maxtime}}};
                result = or_stl_titv(result, temp);
                }     
            }
        }
    }
    for (size_t i = 0; i < sp1.size(); i++)
    {
        for (size_t j = 0; j < sp2.size(); j++)
        {
            if (sp1[i].first.first==1 && sp2[j].first.first==1){
                TimeInterval intersection = computeIntersection(sp1[i], sp2[j], 1); //changer la value
                
                intersection.second.first = intersection.second.first-time_itv.second;
                intersection.second.second = intersection.second.second-time_itv.first;
                intersection = computeIntersection(sp1[i], intersection, 1);
                //intersection.first.second = {-1};
                if (intersection.first.first != -1)
                {
                temp_final = {{{0,{-1}},{0,intersection.second.first}}, intersection, {{0,{-1}}, {intersection.second.second,maxtime}}};
                //print_TimeIntervals(temp_final);
                //std::cout<<"--------"<<std::endl;
                result_final = or_stl_titv(result_final, temp_final);
                }            
            }
            //mettre le backshifting [a,b]
            // mettre l'intersection avec sp1[i]
            //fusionner le resultat dans un vecteur de [0;1] faisant la bonne longueur
        }
    }

    result = or_stl_titv(result_final, result); 
    return result;
}

vector<TimeInterval> Finally(const vector<TimeInterval>& list1, pair<double, double> time_itv){

    vector<TimeInterval> temp = {{{1,{-1}},{0, list1.back().second.second + time_itv.second}}};
    return until_stl_minkowski(temp, list1, time_itv);
}
vector<TimeInterval> Globally(const vector<TimeInterval>& list1, pair<double, double> time_itv){

    return neg_stl_titv(Finally(neg_stl_titv(list1),time_itv));
}

std::vector<std::pair<IntervalVector, Interval>> Sim_to_jn_tube(ibex::simulation& sim){
    std::vector<std::pair<IntervalVector, Interval>> jn_box;
    for (const auto& sol : sim.list_solution_g) {
        if (sol.box_jn && (sol.time_j.lb()>=0)) { // Ensure it's not NULL
            jn_box.push_back({*sol.box_j1,sol.time_j}); // Dereference pointer
        }
    }
    //std::cout<<"conversion ok"<<std::endl;
    return jn_box;
}

std::vector<std::pair<ibex::IntervalVector, ibex::Interval>> Sim_to_jn_tube_test(const ibex::simulation& sim) {
    std::vector<std::pair<ibex::IntervalVector, ibex::Interval>> jn_box;
    jn_box.reserve(sim.list_solution_g.size());  // std::list::size() is valid

    for (auto it = sim.list_solution_g.begin(); it != sim.list_solution_g.end(); ++it) {
        if (it->box_jn && it->time_j.lb() >= 0) {
            jn_box.push_back(std::make_pair(*it->box_j1, it->time_j));
        }
    }

    return jn_box;
}

std::pair<std::vector<vector<TimeInterval>>, std::vector<Interval>> predicate_satisfaction_jn(const std::vector<std::pair<IntervalVector, Interval>>& jn_box, const vector<IntervalVector>& predicate_list) {
    //jn_box;
    int undef_counter = 0;
    std::vector<Interval> undef_titv;
    std::vector<vector<TimeInterval>> result;
    for (size_t j = 0; j < predicate_list.size(); j++)
    {
      std::vector<TimeInterval> P_satisf;
        for (size_t i = 0; i < jn_box.size(); i++)
        {
          int s_val = 0;
          if (jn_box[i].first.is_subset(predicate_list[j]))
          {
            s_val = 1;
          }
          else if (jn_box[i].first.intersects(predicate_list[j]))
          {
            s_val = 2;
          }
          
          if (s_val != 2)
          {
            P_satisf.push_back({{s_val,{-1}},{jn_box[i].second.lb(), jn_box[i].second.ub()}});
          }
          else
          {
            P_satisf.push_back({{s_val,{undef_counter}},{jn_box[i].second.lb(), jn_box[i].second.ub()}});
            undef_counter++;
            undef_titv.push_back(jn_box[i].second); //en connaisant le undef_counter au niveau general on retrouve directement l'interval de temp
          }
      }
      result.push_back(P_satisf);
    }
    //test formule

    //on reboucle
  return {result,undef_titv};
}

std::pair<std::vector<vector<TimeInterval>>, std::vector<Interval>> predicate_satisfaction(ibex::simulation& sim, const vector<IntervalVector>& predicate_list) {
    std::vector<std::pair<IntervalVector, Interval>> jn_box;
    for (const auto& sol : sim.list_solution_g) {
        if (sol.box_jn && (sol.time_j.lb()>=0)) { // Ensure it's not NULL
            jn_box.push_back({*sol.box_j1,sol.time_j}); // Dereference pointer
        }
    }
    //jn_box;
    int undef_counter = 0;
    std::vector<Interval> undef_titv;
    std::vector<vector<TimeInterval>> result;
    for (size_t j = 0; j < predicate_list.size(); j++)
    {
      std::vector<TimeInterval> P_satisf;
        for (size_t i = 0; i < jn_box.size(); i++)
        {
          int s_val = 0;
          if (jn_box[i].first.is_subset(predicate_list[j]))
          {
            s_val = 1;
          }
          else if (jn_box[i].first.intersects(predicate_list[j]))
          {
            s_val = 2;
          }
          
          if (s_val != 2)
          {
            P_satisf.push_back({{s_val,{-1}},{jn_box[i].second.lb(), jn_box[i].second.ub()}});
          }
          else
          {
            P_satisf.push_back({{s_val,{undef_counter}},{jn_box[i].second.lb(), jn_box[i].second.ub()}});
            undef_counter++;
            undef_titv.push_back(jn_box[i].second); //en connaisant le undef_counter au niveau general on retrouve directement l'interval de temp
          }
      }
      result.push_back(P_satisf);
    }
    //test formule

    //on reboucle
  return {result,undef_titv};
}



vector<double> Bisection_time(const vector<Interval>& Un_origin, const double& min_step_time, const double& t_min){
  
   std::vector<double> vec;

  for (size_t i = 0; i < Un_origin.size(); i++)
  {
    if (Un_origin[i].ub()-Un_origin[i].lb() > min_step_time && Un_origin[i].mid()>t_min)
    {
          vec.push_back(Un_origin[i].mid());
    }
  }
    // Trier le vecteur
    std::sort(vec.begin(), vec.end());
    // Supprimer les doublons adjacents
    vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
  return vec;
}

vector<Interval> Uncertainty_origin(vector<TimeInterval> Phi1, const vector<Interval>& P_Satisfaction_signals){
vector<Interval> result;
for (size_t i = 0; i < Phi1.size(); i++)
    {
        if (Phi1[i].first.first == 2)
        {
            Phi1[i].first.second = minusonecleaner(Phi1[i].first.second);
            std::cout << "uncertainty index "<<i<<": ";
            if (Phi1[i].first.second[0]==-2)
            {
              std::cout <<"Incomplete simulation coupling ["<< Phi1[i].second.first <<", "<<Phi1[i].second.second<<"]"<<std::endl;
            }
            
            for (size_t j = 0; j < Phi1[i].first.second.size(); j++)
            {
                if (Phi1[i].first.second[j] >= 0)
                {
                  std::cout << Phi1[i].first.second[j] <<", titv: "<< P_Satisfaction_signals[Phi1[i].first.second[j]] << " | ";
                  result.push_back(P_Satisfaction_signals[Phi1[i].first.second[j]]);
                }
            }
            std::cout<<std::endl;
        }
    }
 return result;
}

int satisfies_at_time(const double& time, const vector<TimeInterval>& phi){
int result;
vector<TimeInterval> Phi1 = phi;
for (size_t i = 0; i < phi.size(); i++)
{
    
  if (phi[i].second.second > time)
  {
    result = phi[i].first.first;
    if (result == 2) ///on vérifie si ça vient de la simu
    {
        Phi1[i].first.second = minusonecleaner(Phi1[i].first.second);
        if (Phi1[i].first.second[0]==-2)
        {
            std::cout <<"Incomplete simulation coupling, must continue ["<< phi[i].second.first <<", "<<phi[i].second.second<<"]"<<std::endl;
            return -2;
        }
    }
    
    return result;

  }  
}
return 3; ///unknown
}

std::vector<std::pair<IntervalVector, Interval>> Append_tube(const std::vector<std::pair<IntervalVector, Interval>>& old_tube, const std::vector<std::pair<IntervalVector, Interval>>& extension_tube){

  std::vector<std::pair<IntervalVector, Interval>>  result = old_tube;
  std::vector<std::pair<IntervalVector, Interval>>  temp_tube = extension_tube;
  double end_time = 0.0;
if (!result.empty()) {
    end_time = result.back().second.ub();
}
             std::cout<<"tube end:"<<end_time<<std::endl;
if (temp_tube.empty())
{
   std::cout<<"c'est la merde"<<std::endl;
}

  for (size_t i = 0; i < temp_tube.size(); i++)
  {
    if (temp_tube[i].second.lb()!=temp_tube[i].second.ub() && temp_tube[i].second.lb()>= 0)
    {   
        Interval temp = Interval(temp_tube[i].second.lb() + end_time, temp_tube[i].second.ub() + end_time);
        temp_tube[i].second = temp;
        result.push_back(temp_tube[i]);
    }
  }
return result;
}

double last_time_tube(const std::vector<std::pair<IntervalVector, Interval>>&  result){
      double end_time = 0.0;
    if (!result.empty()) {
        end_time = result.back().second.ub();
    }
    return end_time;
}

double first_time_tube(const std::vector<std::pair<IntervalVector, Interval>>&  result){
      double end_time = 0.0;
    if (!result.empty()) {
        end_time = result.front().second.lb();

        if (end_time < 0)
        {
            end_time = 0;
        }
        
    }
    else
    {
        std::cerr << "Warning: first_time_tube called with empty tube." << std::endl;
        return 0.0;  // Valeur de secours
    }
    return end_time;
}

void process_and_save_data(ibex::simulation& sim, 
                           const vector<IntervalVector>& predicate_list, 
                           const vector<int>& selected_indices,
                           const string& predicate_file,
                           const string& jn_box_file) {
    // Extract jn_box from sim
    vector<pair<IntervalVector, Interval>> jn_box;
    for (const auto& sol : sim.list_solution_g) {
        if (sol.box_jn && (sol.time_j.lb() >= 0)) { 
            jn_box.push_back({*sol.box_j1, sol.time_j});
        }
    }

    // Open files for writing
    ofstream predFile(predicate_file);
    ofstream jnBoxFile(jn_box_file);

    if (!predFile || !jnBoxFile) {
        cerr << "Error: Unable to open files for writing!" << endl;
        return;
    }

    // Writing selected indices from predicate_list
    for (size_t j = 0; j < predicate_list.size(); j++) {
        for (int n : selected_indices) { // Iterate over selected indices
            if (n < predicate_list[j].size()) { // Ensure index is within bounds
                predFile << predicate_list[j][n] << " ";
            }
        }
        predFile << "\n";
    }
    bool isempty = false;
    // Writing selected indices from jn_box
    for (size_t i = 0; i < jn_box.size(); i++) {
        for (int n : selected_indices) { // Iterate over selected indices
            if (n < jn_box[i].first.size()) { // Ensure index is within bounds
                if (jn_box[i].first[n].is_empty())
                {
                  isempty = true;
                }
                
                if (!isempty)
                {
                  jnBoxFile << jn_box[i].first[n] << " ";
                }
            }
        }
        if (!isempty)
        {
         jnBoxFile << "\n";
        }      
        isempty= false;
    }

    // Close files
    predFile.close();
    jnBoxFile.close();
    std::cout << "Files saved: " << predicate_file << ", " << jn_box_file << endl;
}
