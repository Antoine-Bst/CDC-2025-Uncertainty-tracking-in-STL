/* ============================================================================
 * D Y N I B E X - STL example
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

#include "ibex.h"
#include <iostream>
#include "CSTL.h"
#include <memory>



using namespace ibex;
using namespace std;

struct simu_incr {
  std::vector<std::shared_ptr<ibex::simulation>> Simu_list;
  std::vector<std::pair<IntervalVector, Interval>> jn_tube;
  std::vector<Interval> time_itv;
  std::vector<std::vector<std::pair<IntervalVector, Interval>>> local_jn_tubes;

  // Constructor
  simu_incr(std::shared_ptr<ibex::simulation> sim) {
    Simu_list.push_back(sim);
    jn_tube = Sim_to_jn_tube(*sim);

    double t_start = first_time_tube(jn_tube);
    double t_end = last_time_tube(jn_tube);
    time_itv.push_back(Interval(t_start, t_end));
    local_jn_tubes.push_back(jn_tube);
  }

  // Add new simulation
  void add(std::shared_ptr<ibex::simulation> sim) {
    Simu_list.push_back(sim);

    auto temp_jn_tube = Sim_to_jn_tube(*sim);

    if (temp_jn_tube.empty()) {
      std::cerr << "Warning: empty tube, skipping..." << std::endl;
      time_itv.push_back(time_itv.back());
      local_jn_tubes.push_back({});
      return;
    }

    double last_time = !time_itv.empty() ? time_itv.back().ub() : 0.0;

    double t_start = first_time_tube(temp_jn_tube) + last_time;
    double t_end = last_time_tube(temp_jn_tube) + last_time;

    time_itv.push_back(Interval(t_start, t_end));
    local_jn_tubes.push_back(temp_jn_tube);
  }

  // Contraction
  void contract_at(const double& time) {
    std::cout << "Contracting at time: " << time << std::endl;
    for (size_t i = 0; i < time_itv.size(); ++i) {
      std::cout << " - Segment " << i << " covers [" << time_itv[i].lb() << ", " << time_itv[i].ub() << "]" << std::endl;
    }

    if (time_itv.size() != Simu_list.size()) {
      std::cout << "Size error in incremental simulation object" << std::endl;
      return;
    }

    for (size_t i = 0; i < time_itv.size(); ++i) {
      if (time >= time_itv[i].lb() && time < time_itv[i].ub()) {
        if (!Simu_list[i]) {
          std::cerr << "Null simulation pointer at index " << i << std::endl;
          return;
        }

        double local_time = time - time_itv[i].lb();
        auto before = Sim_to_jn_tube(*Simu_list[i]);

        if (local_time < 0 || local_time > last_time_tube(before)) {
          std::cout << "Invalid contraction time: " << local_time << std::endl;
          return;
        }

        std::cout << "local_bis_time: " << local_time << std::endl;

        Simu_list[i]->get_tight(local_time);

        auto after = Sim_to_jn_tube(*Simu_list[i]);
        std::cout << "Before size: " << before.size() << ", After size: " << after.size() << std::endl;

        local_jn_tubes[i] = after;

        break;
      }
    }
  }

  // Rebuild full tube
  void update_jn_tube() {
    std::vector<std::pair<IntervalVector, Interval>> temp_tube;
    double current_shift = 0.0;

    for (size_t i = 0; i < local_jn_tubes.size(); ++i) {
      const auto& tube = local_jn_tubes[i];

      for (const auto& pair : tube) {
      const IntervalVector& box = pair.first;
      const Interval& t = pair.second;

      if (t.lb() != t.ub() && t.lb() >= 0) {
        Interval shifted_t(t.lb() + current_shift, t.ub() + current_shift);
        temp_tube.emplace_back(box, shifted_t);
      }
    }


      if (!tube.empty()) {
        current_shift += tube.back().second.ub();
      }
    }

    jn_tube = temp_tube;
  }
};



int main(){

  const int n= 4; ///number of state of the differential system
  Variable y(n);
  IntervalVector Box(n);
  IntervalVector yinit(n);
  yinit[0] = Interval(4,4.01);  //x
  yinit[1] = Interval(4,4.01);   //x'
  yinit[2] = Interval(4,4.01);  //x
  yinit[3] = Interval(4,4.01);   //x'
  Affine2Vector yinit_aff(yinit, true);


  vector<IntervalVector> p_list; ///predicate list

  IntervalVector predicate1(n);
    predicate1[0] = Interval(-0.75,0.75);
    predicate1[1] = Interval(-1,1.1);
    predicate1[2] = Interval(-1000,1000);
    predicate1[3] = Interval(-1000,1000);

   IntervalVector predicate2(n);
    predicate2[0] = Interval(1.5,1.9);
    predicate2[1] = Interval(-1.5,-0.8);
    predicate2[2] = Interval(-1000,1000);
    predicate2[3] = Interval(-1000,1000);


    IntervalVector predicate3(n);
    predicate3[0] = Interval(3,4);
    predicate3[1] = Interval(-1.5,0.5);
    predicate3[2] = Interval(-1000,1000);
    predicate3[3] = Interval(-1000,1000);



    p_list = {predicate1,predicate2, predicate3};

    int satisfaction_at_t = -2;
    double timesatisf = 0;

    const double step = 1;
    double tsimu= 1; //initial simulation time (minimum timesatisf)

    Interval K = Interval(1,1.3);
    std::vector<std::pair<IntervalVector, Interval>> result_tube = {};

      Function ydot = Function(y,Return(K*(y[2]-y[0]),K*(y[3]-y[1]),y[3], 0.5*(1-sqr(y[2]))*y[3]-y[2]));
      ivp_ode problem = ivp_ode(ydot,0.0,yinit_aff);    
      auto sim_ptr = std::make_shared<simulation>(&problem, tsimu, KUTTA3, 1e-4, 1e-2); //Time Bisecting
      //auto sim_ptr = std::make_shared<simulation>(&problem, step, KUTTA3, 1e-9, 1e-5); //Only incremental
      sim_ptr->run_simulation();
      yinit_aff = sim_ptr->get_last_aff();
      std::cout<<"simu_init: ok"<<std::endl;
      simu_incr tube_sim(sim_ptr);
      double simulation_time = tsimu;

      while (satisfaction_at_t==-2){ //example with incremental simulation
      //---------------simulation--------------------
      Function ydot = Function(y,Return(K*(y[2]-y[0]),K*(y[3]-y[1]),y[3], 0.5*(1-sqr(y[2]))*y[3]-y[2]));
      ivp_ode problem = ivp_ode(ydot,0.0,yinit_aff);   

      auto new_sim = std::make_shared<simulation>(&problem, step, KUTTA3, 1e-4, 1e-2); //Time Bisecting
      //auto new_sim = std::make_shared<simulation>(&problem, step, KUTTA3, 1e-9, 1e-5); //Only incremental
      
      new_sim->run_simulation();

      std::cout<<"-----starting_loop-----"<<std::endl;
      ///-----------------------fin simu-------------------------
      tube_sim.add(new_sim);
      tube_sim.update_jn_tube();
      
      ////------------------------Verifying STL formula and contracting--------------------------------
        vector<TimeInterval> Phi1;
        vector<double> Bisectime;
        int f = 0;

      
      do
      {
        std::pair<std::vector<vector<TimeInterval>>, std::vector<Interval>> P_Satisfaction_signals = predicate_satisfaction_jn(tube_sim.jn_tube, p_list);
        Phi1 = and_stl_titv(Globally(neg_stl_titv(P_Satisfaction_signals.first[0]), {3, 5}), Globally( or_stl_titv( Finally(P_Satisfaction_signals.first[1], {1,2}), neg_stl_titv(P_Satisfaction_signals.first[2])), {3, 4}));
        //print_TimeIntervals(Phi1);
        std::cout<<"phi1bon"<<std::endl;
        vector<Interval> un_time_itv = Uncertainty_origin(Phi1, P_Satisfaction_signals.second);
        Bisectime = Bisection_time(un_time_itv, 0.05,simulation_time); //the minimum step (0.05 here) should be small if strong and intricated temporal constraints are considered; Put 1 to have only incremental verification

        for (size_t i = 0; i < Bisectime.size(); i++)
        {
          std::cout<< Bisectime[i]<<" , ";
          //std::cout<<"bisect" <<std::endl;
              if (Bisectime[i] < tube_sim.time_itv.front().lb()) {
              std::cout << "Trying to bisect before simulation start — skipped: " << Bisectime[i] << std::endl;
              continue;
              }
          tube_sim.contract_at(Bisectime[i]); //le problème vient qu'il n'a pas le droit de bisecter dans le passé
        }
        f++;
        //std::cout<<"bisect_ok"<<std::endl;
        std::cout<<std::endl;
          auto before = tube_sim.jn_tube;
          tube_sim.update_jn_tube();
          //auto after = tube_sim.jn_tube;
          //std::cout << "Before update size: " << before.size() << ", After update size: " << after.size() << std::endl;
        satisfaction_at_t = satisfies_at_time(timesatisf, Phi1);
        std::cout<<"-------satisfaction à "<<timesatisf<<"s :"<< satisfaction_at_t <<std::endl;
        std::cout<<"----------new_iteration-----------"<< f <<std::endl;
      } while (!Bisectime.empty()&&(satisfaction_at_t==2||satisfaction_at_t==-2));
      

      simulation_time = last_time_tube(tube_sim.jn_tube);

      std::cout<<"----------Current_simulation_time-----------"<< simulation_time <<std::endl;
      //////----------------------end adaptative----------------------------
      yinit_aff = new_sim->get_last_aff();
      }
      
  return 0;
}

