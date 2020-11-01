**AMORE - DSI Capstone FALL 2020**
**FILES**
- **THIS README FILE**
- **main.py**
  - read in data
  - calculate weights
- **isoprene_rates.py** - functions used to calculate weights
  - TUN, ALK, NIT, ESO1, ESO2, TROE, FALL (adding more functions if necessary, see the rtf file)
  - Other supplemental math functions
- **read_input.py**   - read in spc, def, eqn files
  - read_spc(file)
  - read_eqn(file)
  - read_def(file)
- **calculations.py** 
  - calculate_weight(eqn, init_values): return weight dict  
  - **TODO: calculate_all_weights(eqns, init_values): return a weight dict**
  - **TODO: calculate_r**
  - **TODO: calculate_all_r: return a r_ab dict**
- **directed_graph.py**  
  - **TODO: DFS
  - **TODO: find_dependent_set
- **requirements.txt**
  - pip install -r requirements.txt

**TODO:**  
*calculations.py*  
Direct influence - a normalized contribution of species B to the production rate of species A, rAB, if the normalized contribution r_AB is sufficiently large, species A strongly depends on species B.  
arguments might need to be changed in this function, because we have the a weight_dict from calculate_weight function  
Define calculate_r(v_a vector, w vector, d_b vector):  

CLASS 2 - *directed_graph.py*  
FUNCTION 1  
Define get_nodes(all r_AB, threshold epsilon):    
Use the threshold epsilon and compare it with each of the r_AB (direct influence) to get a subset of nodes that are valid and important.  
Return (node set (A, B, etc))  

FUNCTION 2  
Define get_edges(all r_AB, threshold epsilon)    
Use the threshold epsilon and compare it with each of the r_AB (direct influence) to get a subset of edges for DRG construction.  
Return (edge set (A -> B))  

FUNCTION 3
after getting nodes and edges  
Define construct_DRG (nodes, edges, threshold epsilon):  
the DRG can be constructed by the following rules:  
(1) Each node in DRG is uniquely mapped to a species in the detailed mechanism.  
(2) There exists a directed edge from A to B if and only if r_AB is larger than or equal to epsilon.  

FUNCTION 4  
Define find_dependent_set(self/graph, species A):  
DFS  

FUNCTION 5 (hold)  
insert species into DRG  
Define insert_species(graph/self, species A, threshold epsilon):  
Species A is inserted into the DRG by adding a note to represent the species. Then based on the threshold, we use the threshold to add edges.  


