Welcome ! This is a prototype for the verification of STL formula on reachable tube using Dynibex. It has been tested on Linux Ubuntu only.
It uses DynIbex.

Install Dynibex:
https://perso.ensta-paris.fr/~chapoutot/dynibex/index.php#download-installation

In recent Ubuntu versions you might need to do the the install in a python 2.7 virtual environement and :
sudo CXXFLAGS="-std=c++14" ./waf configure
sudo CXXFLAGS="-std=c++14" ./waf install

If dynibex is a local install add to the make file:
export PKG_CONFIG_PATH='path_to_dynibex'/share/pkgconfig 

Make the code:
Open terminal in the CSTL directory and compile using make:
make

Run:
./simulation.out

If everything works correctly, the output should be for the last two line:
-------satisfaction à 0s :1
----------new_iteration-----------1
----------Current_simulation_time-----------6


=====================================Small Documentation===================================
The verification is performed in a bottom-up manner using the syntax tree of the formula and the satisfaction signals

phi1 = neg_stl(phi);             // Logical negation: ¬phi
phi1 = and_stl(phi2, phi3);      // Logical AND: phi2 ∧ phi3
phi1 = or_stl(phi2, phi3);       // Logical OR: phi2 ∨ phi3
phi1 = until_stl(phi2, phi3, {t1, t2});  // Until operator: phi2 U[t1,t2] phi3
phi1 = Finally(phi, {t1, t2});   // Eventually operator: F[t1,t2] phi
phi1 = Globally(phi, {t1, t2});  // Always operator: G[t1,t2] phi


predicate_satisfaction(sim, predicates); Constructs satisfaction signals for the list of predicates and the simulation object.
The output is a list of signals corresponding to each predicate in the input list.

print_Satisf_Signals(phi); Displays the satisfaction signal of a given formula.
