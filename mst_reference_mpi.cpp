#include <vector>
#include <mpi.h>
#include <float.h>
#include <stdlib.h>
#include <assert.h>
#include "defs.h"
#include <algorithm>
// #include <cstdint>
#include <cstring>
#include <iostream>
#include <omp.h>
#include <unistd.h>

#include <map>

#define UINT32_MAX  (0xffffffff)
#define UINT64_MAX  (0xffffffffffffffff)



using namespace std;
// map <vertex_id_t,vertex_id_t> components_map;
typedef map <vertex_id_t, loc_edge_t>  mapT;

void print_source_graph(graph_t *G)
{
  for(int p=0; p < G->nproc; p++){
    if(p == G->rank){
      printf("\nnproc=%d n=%d m=%d\n", G->nproc, G->local_n, G->local_m);
      printf("RANK= %d \n", G->rank);
      for(vertex_id_t i = 0; i < G->local_n; i++) {
        printf("i= %d|||| \n", VERTEX_TO_GLOBAL(i, G->n, G->nproc, G->rank));
        for (vertex_id_t j = G->rowsIndices[i]; j < G->rowsIndices[i+1]; j++) {
          int a = VERTEX_TO_GLOBAL(i, G->n, G->nproc, G->rank);
            printf("%d to %d, weight= %lf, edge_number= %d  \n", a, G->endV[j], G-> weights[j], edge_to_global(j,G));
        }

      }
      printf("\n");
    }
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

/* returns global number of edge */
edge_id_t edge_to_global(edge_id_t edge, graph_t *G) {
    int rank = G->rank;
    int size = G->nproc;
    edge_id_t g_edge = 0;
    for(int i = 0; i < rank && i < size; ++i)
    {
        g_edge += G->num_edges_of_any_process[i];
    }
    return (g_edge + edge);
}

/* NOTE: isolated vertex is also tree, such tree must be represented as separate element of trees_mst vector with zero-length edges list
 * FIXME: If you change MST output data structure, you must change this function */
typedef vector<vector<edge_id_t > > result_t;

extern "C" void convert_to_output(graph_t *G, result_t &trees_mst, forest_t *trees_output)
{
    int MPI_TAG = 99;
    int number_of_edges = 0;
    map <vertex_id_t,vertex_id_t> components_map;
    // result_t &trees_mst = *reinterpret_cast<result_t*>(result);
    vertex_id_t trees_count = trees_mst.size();
    edge_id_t max_number = 0;
    for(vertex_id_t i = 0; i < G->nproc; i++){
      if(max_number < G->num_edges_of_any_process[i]){
        max_number = G->num_edges_of_any_process[i];
      }
    }

    edge_id_t* buf_edges = new edge_id_t[max_number];
    if(G->rank == 0){
      vertex_id_t components_count = 0;
      for(vertex_id_t i = 0; i < trees_count; i++)
      {
        vertex_id_t component =  trees_mst[i][0];
        number_of_edges += trees_mst[i].size() - 1;
        map <vertex_id_t,vertex_id_t>::iterator found = components_map.find(component);
        if(found == components_map.end()){
          components_map.insert(pair<vertex_id_t, vertex_id_t>(component, components_count));
          // printf("component= %u\n", component);
          components_count++;
        }
      }


      for(vertex_id_t i = 0; i < G->nproc-1; i++){
        vertex_id_t trees_size;
        MPI_Status status;
        MPI_Recv(&trees_size, 1,MPI_UINT32_T, MPI_ANY_SOURCE, MPI_TAG, MPI_COMM_WORLD, &status);
        int proc = status.MPI_SOURCE;
        for(vertex_id_t j = 0; j < trees_size; j++){
          MPI_Recv(buf_edges, max_number,MPI_UINT64_T, proc, MPI_TAG, MPI_COMM_WORLD, &status);
          vertex_id_t component = buf_edges[0];
          int count;
          MPI_Get_count(&status, MPI_INT, &count);
          count = count / 2;
          number_of_edges+= count - 1;
          map <vertex_id_t,vertex_id_t>::iterator found = components_map.find(component);
          if(found == components_map.end()){
            components_map.insert(pair<vertex_id_t, vertex_id_t>(component, components_count));
            trees_mst.push_back(vector<edge_id_t>());
            trees_mst[components_count].push_back(component);
            for(int k = 1; k < count; k++)
            {
              trees_mst[components_count].push_back(buf_edges[k]);
            }
            components_count++;
          }else{
            for(int k = 1; k < count; k++)
            {
              trees_mst[found->second].push_back(buf_edges[k]);
            }

          }
        }
      }
      trees_output->p_edge_list = (edge_id_t *)malloc(trees_mst.size()*2 * sizeof(edge_id_t));
      trees_output->edge_id = (edge_id_t *)malloc(number_of_edges * sizeof(edge_id_t));
      trees_output->p_edge_list[0] = 0;
      trees_output->p_edge_list[1] = trees_mst[0].size() - 1;
      for (vertex_id_t i = 1; i < trees_mst.size(); i++) {
          trees_output->p_edge_list[2*i] = trees_output->p_edge_list[2*i-1];
          trees_output->p_edge_list[2*i +1] = trees_output->p_edge_list[2*i-1] + trees_mst[i].size() - 1;
      }
      int k = 0;
      for (vertex_id_t i = 0; i < trees_mst.size(); i++) {
          for (edge_id_t j = 1; j < trees_mst[i].size(); j++) {
              trees_output->edge_id[k] = trees_mst[i][j];
              k++;
          }
      }
      printf("size= %d\n",  G->nproc);
      printf("numTrees= %u\n",  trees_mst.size());
      fflush(stdout);
      trees_output->numTrees = trees_mst.size();
      trees_output->numEdges = number_of_edges;

    }else{
      MPI_Send(&trees_count, 1,MPI_UINT32_T,0, MPI_TAG, MPI_COMM_WORLD);
      for(vertex_id_t i = 0; i < trees_count; i++)
      {
          MPI_Send(&trees_mst[i][0], trees_mst[i].size(),MPI_UINT64_T, 0, MPI_TAG, MPI_COMM_WORLD);

      }
    }
    free(buf_edges);

}

extern "C" void init_mst(graph_t *G) {
    G->num_edges_of_any_process = (edge_id_t*) malloc (G->nproc * sizeof(edge_id_t));
    edge_id_t * edges_to_send = (edge_id_t*) malloc (G->nproc * sizeof(edge_id_t));
    assert(G->num_edges_of_any_process);

    for(int i = 0; i < G->nproc; ++i)
  edges_to_send[i] = G->local_m;

    MPI_Alltoall(edges_to_send, 1, MPI_UINT64_T, G->num_edges_of_any_process, 1, MPI_UINT64_T, MPI_COMM_WORLD);
    if(edges_to_send) free(edges_to_send);
}

vertex_id_t find(full_edge_t* edges, vertex_id_t i, graph_t* G)
{
  vertex_id_t new_i;
  while(i != edges[i].parent){
    new_i = edges[i].parent;
    edges[i].parent = edges[new_i].parent;
    i = edges[new_i].parent;
  }
  return i;
}

edge_id_t find_in_edges(edge_id_t* parents, edge_id_t i, graph_t* G)
{
  edge_id_t new_i;
  while(parents[i] != i){
    new_i = parents[i];
    parents[i] = parents[new_i];
    i = parents[new_i];
  }
  return i;
}

vertex_id_t component_from_find(full_edge_t* edges, vertex_id_t i, vertex_id_t first_vertex,
                                graph_t* G, edge_id_t* edge_parents, map <vertex_id_t, edge_id_t> &edges_map){
  vertex_id_t from_find = find(edges, i, G);
  vertex_id_t foreign = edges[from_find].foreign;
  vertex_id_t component;
  if(foreign == UINT32_MAX){
    component = first_vertex + from_find;
  }else{
    component = G->endV[find_in_edges(edge_parents,edges_map.find(foreign)->second,G)];
    edges[i].foreign = component;
  }
  return component;
}


void main_merge(full_edge_t* full_edges,vertex_id_t* broadcast_results, graph_t* G,
               edge_id_t *edge_parents, map <vertex_id_t, loc_edge_t> &edges,
               bool &changed, bool *added_edges){
  vertex_id_t first_vertex = VERTEX_TO_GLOBAL(0, G->n, G->nproc, G->rank);
  vertex_id_t results_count = 0;
  for(vertex_id_t i = 0;i < G->local_n;i++){
    vertex_id_t from = i + first_vertex;
    if(full_edges[i].foreign != from)
      continue;

    vertex_id_t to = full_edges[i].to;
    vertex_id_t cheapest_received_from = full_edges[i].cheapest_rank_to;
    mapT::iterator found_to = edges.find(to);
    mapT::iterator found_from = edges.find(from);
    if(found_to != edges.end() && found_to->second.edge != UINT64_MAX ){
      edge_id_t edge = find_in_edges(edge_parents, found_to->second.edge, G);
      if(G->endV[edge] != to){
        to = G->endV[edge];
        found_to = edges.find(to);
      };
    }

    if(to == UINT32_MAX){
      continue;
    }
    if(from == to){
      broadcast_results[results_count] = from;
      broadcast_results[results_count + 1] = to;
      broadcast_results[results_count + 2] = cheapest_received_from;
      results_count+=3;
      continue;
    }

    if(found_from == edges.end()){
      loc_edge_t loc_edge;
      loc_edge.local = UINT32_MAX;
      loc_edge.edge = UINT64_MAX;
      if(found_to != edges.end()){
        edges.insert( pair<vertex_id_t, loc_edge_t>(from, loc_edge) );
        found_from = edges.find(from);
      }else{
        continue;
      }
    }
    if(found_to == edges.end()){
      loc_edge_t loc_edge;
      loc_edge.local = found_from->second.local;
      loc_edge.edge = found_from->second.edge;

      edges.insert( pair<vertex_id_t, loc_edge_t>(to, loc_edge) );
      if(loc_edge.edge != UINT64_MAX)
        G->endV[loc_edge.edge] = to;
      if(loc_edge.local != UINT32_MAX){
        full_edges[loc_edge.local].foreign = to;
        found_from->second.local = UINT32_MAX;
      }

    }else{
      vertex_id_t from_local = found_from->second.local;
      vertex_id_t to_local = found_to->second.local;
      edge_id_t from_edge = found_from->second.edge;
      edge_id_t to_edge = found_to->second.edge;
      if(from_local != UINT32_MAX){
        if(G->rank == cheapest_received_from){
          added_edges[full_edges[from_local].edge] = true;
        }
        if(to_local == UINT32_MAX){
          full_edges[from_local].foreign = to;
          found_to->second.local = from_local;

        }else{
          full_edges[from_local].parent = to_local;
          full_edges[from_local].foreign = UINT32_MAX;
          if(full_edges[from_local].weight < full_edges[to_local].weight){
            full_edges[to_local].to = full_edges[from_local].to;
            full_edges[to_local].weight = full_edges[from_local].weight;
            full_edges[to_local].edge = full_edges[from_local].edge;
            full_edges[to_local].cheapest_rank_to = full_edges[from_local].cheapest_rank_to;
          }
        }
        found_from->second.local = UINT32_MAX;
      }

      if(from_edge == UINT64_MAX){
        if(to_edge == UINT64_MAX){
        }else{
          found_from->second.edge = to_edge;
        }
      }else{
        if(to_edge == UINT64_MAX){
          G->endV[from_edge] = to;
          found_to->second.edge = from_edge;

        }else{
          G->endV[from_edge] = to;
          edge_parents[from_edge] = found_to->second.edge;
        }
      }
    }


    broadcast_results[results_count] = from;
    broadcast_results[results_count + 1] = to;
    broadcast_results[results_count + 2] = cheapest_received_from;
    results_count+=3;


  }
  MPI_Bcast(&results_count, 1, MPI_UINT32_T, G->rank, MPI_COMM_WORLD);
  MPI_Bcast(broadcast_results, results_count, MPI_UINT32_T, G->rank, MPI_COMM_WORLD);
  if(results_count > 0){
    changed = true;
  }
}





void merge_ask(full_edge_t* full_edges,vertex_id_t* broadcast_results, vertex_id_t results_count, graph_t* G,
               edge_id_t *edge_parents, map <vertex_id_t, loc_edge_t> &edges,
               bool &changed, bool *added_edges){
  vertex_id_t first_vertex = VERTEX_TO_GLOBAL(0, G->n, G->nproc, G->rank);


  if(results_count > 0){
    changed = true;
  }

  for(vertex_id_t i = 0;i < results_count;i+=3){
    vertex_id_t from = broadcast_results[i];
    vertex_id_t to = broadcast_results[i+1];


    vertex_id_t cheapest_received_from = broadcast_results[i+2];

    mapT::iterator found_to = edges.find(to);

    mapT::iterator found_from = edges.find(from);
    if(to == UINT32_MAX){
      continue;
    }
    if(from == to){
      continue;
    }
    if(found_from == edges.end()){
      continue;
    }
    if(found_to == edges.end()){
      loc_edge_t loc_edge;
      loc_edge.local = found_from->second.local;
      loc_edge.edge = found_from->second.edge;

      edges.insert( pair<vertex_id_t, loc_edge_t>(to, loc_edge) );
      if(loc_edge.edge != UINT64_MAX)
        G->endV[loc_edge.edge] = to;
      if(loc_edge.local != UINT32_MAX){
        full_edges[loc_edge.local].foreign = to;
        found_from->second.local = UINT32_MAX;
      }

    }else{
      vertex_id_t from_local = found_from->second.local;
      vertex_id_t to_local = found_to->second.local;
      edge_id_t from_edge = found_from->second.edge;
      edge_id_t to_edge = found_to->second.edge;
      if(from_local != UINT32_MAX){
        if(G->rank == cheapest_received_from){
          added_edges[full_edges[from_local].edge] = true;
        }
        if(to_local == UINT32_MAX){
          full_edges[from_local].foreign = to;
          found_to->second.local = from_local;

        }else{
          full_edges[from_local].parent = to_local;
          full_edges[from_local].foreign = UINT32_MAX;
          if(full_edges[from_local].weight < full_edges[to_local].weight){
            full_edges[to_local].to = full_edges[from_local].to;
            full_edges[to_local].weight = full_edges[from_local].weight;
            full_edges[to_local].edge = full_edges[from_local].edge;
            full_edges[to_local].cheapest_rank_to = full_edges[from_local].cheapest_rank_to;
          }
        }
        found_from->second.local = UINT32_MAX;
      }

      if(from_edge == UINT64_MAX){
        if(to_edge == UINT64_MAX){
        }else{
          found_from->second.edge = to_edge;
        }
      }else{
        if(to_edge == UINT64_MAX){
          G->endV[from_edge] = to;
          found_to->second.edge = from_edge;

        }else{
          G->endV[from_edge] = to;
          edge_parents[from_edge] = found_to->second.edge;
        }
      }
    }
  }
}


extern "C" void MST_boruvka(graph_t *G, result_t &trees) {
    vertex_id_t max_vertex_number = G->n / G->nproc + 1;
    to_weight_t *to_weights = new to_weight_t[max_vertex_number];
    vertex_id_t *broadcast_results = new vertex_id_t[3*max_vertex_number];
    int rank = G->rank, size = G->nproc;
    full_edge_t *full_edges = new full_edge_t[G->local_n];
    edge_id_t *edge_parents = new edge_id_t[G->local_m];
    map <vertex_id_t, loc_edge_t> edges;
    bool *added_edges = new bool[G->local_m];
    memset(added_edges,0,G->local_m * sizeof(bool));

    int weight_nitems = 3;
    int weight_blocklengths[] = {1,1,1};
    MPI_Datatype weight_types[] = {MPI_UINT32_T, MPI_UINT32_T, MPI_DOUBLE};
    MPI_Datatype mpi_weight;
    MPI_Aint     weight_offsets[weight_nitems];

    weight_offsets[0] = offsetof(to_weight_t, from_component);
    weight_offsets[1] = offsetof(to_weight_t, to);
    weight_offsets[2] = offsetof(to_weight_t, weight);

    MPI_Type_create_struct(weight_nitems, weight_blocklengths, weight_offsets, weight_types, &mpi_weight);
    MPI_Type_commit(&mpi_weight);
    vertex_id_t first_vertex = VERTEX_TO_GLOBAL(0, G->n, G->nproc, G->rank);
    for(edge_id_t i = 0; i < G->local_m; i++){
      vertex_id_t vertex = G->endV[i];
      mapT::iterator found = edges.find(vertex);
      if(found == edges.end()){
        loc_edge_t loc_edge;
        loc_edge.edge = i;
        loc_edge.local = UINT32_MAX;

        edges.insert( pair<vertex_id_t,loc_edge_t>(vertex,loc_edge) );
        edge_parents[i] = i;
      }else{
        edge_parents[i] = found->second.edge;
      }
    }
    for(vertex_id_t i = 0; i < G->local_n; i++){
      full_edges[i].parent = i;
      full_edges[i].foreign = first_vertex + i;
      mapT::iterator found = edges.find(i+first_vertex);
      if(found == edges.end()){
        loc_edge_t loc_edge;
        loc_edge.local = i;
        loc_edge.edge = UINT64_MAX;

        edges.insert( pair<vertex_id_t,loc_edge_t>(i+first_vertex,loc_edge) );
      }else{
        found->second.local = i;
      }
    }
    int iter = 0;
    while(true){
        MPI_Status status;
        bool changed = false;
        for(vertex_id_t i = 0; i < G->local_n; i++){
          full_edges[i].to = UINT32_MAX;
          full_edges[i].weight = DBL_MAX;
          full_edges[i].cheapest_rank_to = G->rank;
        }
        for(vertex_id_t i = 0; i < G->local_n; i++) {
          vertex_id_t first_component_index =  find(full_edges, i, G);
          vertex_id_t global_first_component_index = full_edges[first_component_index].foreign;
          for (edge_id_t j = G->rowsIndices[i]; j < G->rowsIndices[i+1]; j++) {
            vertex_id_t second_component_index = G->endV[find_in_edges(edge_parents,j,G)];

            if(global_first_component_index == second_component_index){
              continue;
            }
            if(full_edges[first_component_index].weight > G->weights[j]){
              full_edges[first_component_index].to = second_component_index;
              full_edges[first_component_index].weight = G->weights[j];
              full_edges[first_component_index].edge = j;
            }
          }
        }
        mapT::iterator it = edges.begin();
        for(int p = 0; p < G->nproc; p++) {
          if(p == G->rank){
            for(int p2 = 0; p2 < G->nproc - 1; p2++) {
              MPI_Recv(to_weights, max_vertex_number,mpi_weight, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
              vertex_id_t cheapest_rank_to = status.MPI_SOURCE;
              int count;
              MPI_Get_count(&status, MPI_INT, &count);
              count = count / 4;
              for(vertex_id_t i = 0; i < count; i++) {
                to_weight_t * to_weight =  &to_weights[i];
                vertex_id_t local_from = to_weight->from_component - first_vertex;
                if(full_edges[local_from].weight > to_weight -> weight){
                  vertex_id_t to = to_weight -> to;
                  full_edges[local_from].weight = to_weight -> weight;
                  full_edges[local_from].to = to;
                  full_edges[local_from].cheapest_rank_to = cheapest_rank_to;
                }
              }
            }
            vertex_id_t count = 0;
            main_merge(full_edges, broadcast_results, G, edge_parents, edges, changed, added_edges);
          }
          else{

            vertex_id_t count = 0, results_count;
            for(;VERTEX_OWNER(it->first, G->n, G->nproc) < p && it != edges.end();it++){
            }
            if(VERTEX_OWNER(it->first, G->n, G->nproc) != p || it == edges.end()){
              MPI_Send(to_weights, 0,mpi_weight,p, 0,MPI_COMM_WORLD);
            }else{
              for(; it != edges.end() && VERTEX_OWNER(it->first, G->n, G->nproc) == p ;it++){
                if(it->second.local != UINT32_MAX && full_edges[it->second.local].weight != DBL_MAX){

                  vertex_id_t local_from = it->second.local;
                  vertex_id_t to = full_edges[local_from].to;
                  mapT::iterator found = edges.find(to);
                  if(found->second.edge != UINT64_MAX)
                    to = G->endV[find_in_edges(edge_parents,found->second.edge,G)];
                  to_weights[count].to = to;
                  to_weights[count].from_component = full_edges[local_from].foreign;
                  to_weights[count].weight = full_edges[local_from].weight;

                  count++;
                }
              }
              MPI_Send(to_weights, count,mpi_weight,p, 0,MPI_COMM_WORLD);
            }

            MPI_Bcast(&results_count, 1, MPI_UINT32_T, p, MPI_COMM_WORLD);
            MPI_Bcast(broadcast_results, results_count, MPI_UINT32_T, p, MPI_COMM_WORLD);
            merge_ask(full_edges, broadcast_results, results_count, G, edge_parents, edges, changed, added_edges);
          }
        }
    if(!changed){
      map <vertex_id_t,vertex_id_t> components_map;
      trees.clear();
      vertex_id_t count = 0;
      for(vertex_id_t i = 0; i < G->local_n;i++){
        vertex_id_t global_vertex = i + first_vertex;
        vertex_id_t component_index =  find(full_edges, i, G);
        vertex_id_t foreign_index =  full_edges[component_index].foreign;
        map <vertex_id_t,vertex_id_t>::iterator found = components_map.find(foreign_index);
        vertex_id_t address = 0;
        if(found == components_map.end()){
          components_map.insert(pair<vertex_id_t, vertex_id_t>(foreign_index, count));
          trees.push_back(vector<edge_id_t>());
          trees[count].push_back(foreign_index);
          address = count;
          count++;
        }else{
          address = found->second;
        }



        for(edge_id_t j = G->rowsIndices[i]; j < G->rowsIndices[i+1];j++){
          if(added_edges[j]){
            trees[address].push_back(edge_to_global(j,G));
          }
        }
      }
      free(full_edges);
      free(added_edges);
      free(edge_parents);
      free(to_weights);
      free(broadcast_results);
      return;
    }
  }
}


extern "C" void* MST(graph_t *G) {
    result_t trees;

    trees.clear();
    int rank = G->rank, size = G->nproc;
    vertex_id_t TotVertices = G->n;

    // marked edges are those that lead to vertices already in the tree
    vector<uint8_t> marked_edges(G->local_m, 0);
    // marked vertices are local edges already in the tree
    vector<uint8_t> marked_vertices(G->local_n, 0);

    // start with first vertex on first node
    vertex_id_t root = 0;
    do {
        // start a new tree
        trees.push_back(vector<edge_id_t>());
        // local queue of vertices
        vector<vertex_id_t> queue;

        // keep track of last added edge to mark edges
        vertex_id_t last_vertex = root;

        while (true) {
            if (VERTEX_OWNER(last_vertex, TotVertices, size) == G->rank) {
                // last vertex is ours - put it in the queue and mark
                vertex_id_t last_local_vertex = VERTEX_LOCAL(last_vertex, TotVertices, size, rank);
                marked_vertices[last_local_vertex] = 1;
                queue.push_back(last_local_vertex);
            }

            // mark edges that lead to the last added vertex
            for (edge_id_t j = 0; j < G->local_m; j++) {
                if (G->endV[j] == last_vertex) {
                    marked_edges[j] = 1;
                }
            }

            // determine our best candidate edge
            struct {
                double weight;
                int rank;
                edge_id_t edge;
            } best;
            best.weight = DBL_MAX;
            for (vertex_id_t i = 0; i < queue.size(); i++) {
                for (edge_id_t j = G->rowsIndices[queue[i]]; j < G->rowsIndices[queue[i] + 1]; j++) {
                    // skip marked edges
                    if (!marked_edges[j]) {
                        // check if this edge is better than what we have up to now
                        if (best.weight == DBL_MAX || G->weights[j] < best.weight) {
                            best.weight = G->weights[j];
                            best.edge = j;
                        }
                    }
                }
            }

            // reduce and determine global best edge
            best.rank = G->rank;
            MPI_Allreduce(MPI_IN_PLACE, &best, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
            if (best.weight == DBL_MAX) {
                // no suitable edge found, finish this tree
                break;
            } else {
                if (best.rank == G->rank) {
                    // we have the best edge
                    trees.back().push_back(edge_to_global(best.edge,G));
                    last_vertex = G->endV[best.edge];
                }
                MPI_Bcast(&last_vertex, 1, MPI_UINT32_T, best.rank, MPI_COMM_WORLD);
            }
        }

        // find root of a new tree
        root = UINT32_MAX;
        for (vertex_id_t i = 0; i < G->local_n; i++) {
            if (!marked_vertices[i]) {
                root = VERTEX_TO_GLOBAL(i, TotVertices, size, rank);
                break;
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, &root, 1, MPI_UINT32_T, MPI_MIN, MPI_COMM_WORLD);
    } while (root!=UINT32_MAX);
    return &trees;
}

extern "C" void finalize_mst(graph_t* G) {
   if (G->num_edges_of_any_process) free(G->num_edges_of_any_process);
}
