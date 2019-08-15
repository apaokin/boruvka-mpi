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
#define MPI_UINT8_T  MPI_CHAR
#define MPI_C_BOOL  MPI_CHAR
#define MPI_UINT32_T MPI_UNSIGNED
#define MPI_UINT64_T MPI_UNSIGNED_LONG



using namespace std;
// map <vertex_id_t,vertex_id_t> components_map;

//FIRST
//
//v1 want to connect to v2
// v1
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
result_t trees;

extern "C" void convert_to_output(graph_t *G, void *result, forest_t *trees_output)
{
    int MPI_TAG = 99;
    int number_of_edges = 0;
    map <vertex_id_t,vertex_id_t> components_map;
    result_t &trees_mst = *reinterpret_cast<result_t*>(result);
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
        auto found = components_map.find(component);
        if(found == components_map.end()){
          components_map.insert(pair<vertex_id_t, vertex_id_t>(component, components_count));
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
          auto found = components_map.find(component);
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

      trees_output->numTrees = trees_mst.size();
      trees_output->numEdges = number_of_edges;
      free(buf_edges);

    }else{
      MPI_Send(&trees_count, 1,MPI_UINT32_T,0, MPI_TAG, MPI_COMM_WORLD);
      for(vertex_id_t i = 0; i < trees_count; i++)
      {
          MPI_Send(&trees_mst[i][0], trees_mst[i].size(),MPI_UINT64_T, 0, MPI_TAG, MPI_COMM_WORLD);

      }
    }
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

vertex_id_t merge_ask(full_edge_t* full_edges, vertex_id_t from, graph_t* G, int p,
               edge_id_t *edge_parents, map <vertex_id_t, loc_edge_t> &edges,
               bool &changed, vertex_id_t cheapest_received_from = 0){
  vertex_id_t first_vertex = VERTEX_TO_GLOBAL(0, G->n, G->nproc, G->rank);
  vertex_id_t to;
  vertex_id_t to_with_source[2];
  if(p == G->rank){
    vertex_id_t local_from = from - first_vertex;
    to = full_edges[local_from].to;
    auto found = edges.find(to);
    if(found != edges.end() || found->second.edge != UINT64_MAX ){
      edge_id_t edge = find_in_edges(edge_parents, found->second.edge, G);
      to = G->endV[edge];
      full_edges[local_from].to = to;

    }
    if(full_edges[local_from].weight == DBL_MAX){
      to = UINT32_MAX;
    }


    to_with_source[0] = to;
    to_with_source[1] = cheapest_received_from;

    MPI_Bcast( to_with_source, 2, MPI_UINT32_T, p, MPI_COMM_WORLD );
  }
  else{
    MPI_Bcast( to_with_source, 2, MPI_UINT32_T, p, MPI_COMM_WORLD );
    to = to_with_source[0];
  }
  if(from == to || to == UINT32_MAX){
    return UINT32_MAX;
  }
  changed = true;
  auto found_from = edges.find(from);

  if(found_from == edges.end()){
    return UINT32_MAX;
  }
  auto found_to = edges.find(to);

  if(found_to == edges.end()){
    loc_edge_t loc_edge;
    loc_edge.local = found_from->second.local;
    loc_edge.edge = found_from->second.edge;

    edges.insert( pair<vertex_id_t, loc_edge_t>(to, loc_edge) );
    if(loc_edge.edge != UINT64_MAX)
      G->endV[loc_edge.edge] = to;
    if(loc_edge.local != UINT32_MAX)
      full_edges[loc_edge.local].foreign = to;

  }else{
    vertex_id_t from_local = found_from->second.local;
    vertex_id_t to_local = found_to->second.local;
    edge_id_t from_edge = found_from->second.edge;
    edge_id_t to_edge = found_to->second.edge;

    if(from_local != UINT32_MAX){
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
        }
      }
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
        found_to->second.edge = to_edge;
        edge_parents[from_edge] = found_to->second.edge;

      }
    }
  }
  return to_with_source[1];
}


extern "C" void* MST_boruvka(graph_t *G) {
    // printf("MPI-%d\n", MPI_UINT32_T);
    int size_type;
    // if(G->rank == 0){
    //   MPI_Type_size(MPI_CHAR, &size_type);
    //   printf("%d\n", size_type);
    //   MPI_Type_size(MPI_UNSIGNED, &size_type);
    //   printf("%d\n", size_type);
    //   MPI_Type_size(MPI_SHORT, &size_type);
    //   printf("%d\n", size_type);
    //   MPI_Type_size(MPI_UNSIGNED_LONG, &size_type);
    //   printf("%d\n", size_type);
    //
    // }

    int rank = G->rank, size = G->nproc;
    bool changed;
    full_edge_t *full_edges = new full_edge_t[G->local_n];
    edge_id_t *edge_parents = new edge_id_t[G->local_m];
    map <vertex_id_t, loc_edge_t> edges;
    bool *added_edges = new bool[G->local_m];
    memset(added_edges,0,G->local_m * sizeof(bool));

    int weight_nitems = 2;
    int weight_blocklengths[] = {1,1};
    MPI_Datatype weight_types[] = {MPI_UINT32_T, MPI_DOUBLE};
    MPI_Datatype mpi_weight;
    MPI_Aint     weight_offsets[weight_nitems];

    weight_offsets[0] = offsetof(to_weight_t, to);
    weight_offsets[1] = offsetof(to_weight_t, weight);

    MPI_Type_create_struct(weight_nitems, weight_blocklengths, weight_offsets, weight_types, &mpi_weight);
    MPI_Type_commit(&mpi_weight);


    vertex_id_t first_vertex = VERTEX_TO_GLOBAL(0, G->n, G->nproc, G->rank);
    for(edge_id_t i = 0; i < G->local_m; i++){
      vertex_id_t vertex = G->endV[i];
      auto found = edges.find(vertex);
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
      auto found = edges.find(i+first_vertex);
      if(found == edges.end()){
        loc_edge_t loc_edge;
        loc_edge.local = i;
        loc_edge.edge = UINT64_MAX;

        edges.insert( pair<vertex_id_t,loc_edge_t>(i+first_vertex,loc_edge) );
      }else{
        found->second.local = i;
      }
    }
    while(true){
        bool changed = false;
        for(vertex_id_t i = 0; i < G->local_n; i++){
          full_edges[i].to = UINT32_MAX;
          full_edges[i].weight = DBL_MAX;

        }
        for(vertex_id_t i = 0; i < G->local_n; i++) {
          for (edge_id_t j = G->rowsIndices[i]; j < G->rowsIndices[i+1]; j++) {
            vertex_id_t first_component_index =  find(full_edges, i, G);
            vertex_id_t global_first_component_index = full_edges[first_component_index].foreign;
            vertex_id_t second_component_index = G->endV[find_in_edges(edge_parents,j,G)];
            if(first_component_index == UINT32_MAX){
            }

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
        for(int p = 0; p < G->nproc; p++) {
          if(p == G->rank){
            for(vertex_id_t i = 0; i < G->local_n; i++) {
              full_edge_t *cur_full_edge = &full_edges[i];
              if(cur_full_edge->foreign !=UINT32_MAX && cur_full_edge->foreign == i + first_vertex){
                vertex_id_t cur_vertex = first_vertex + i;
                MPI_Bcast( &cur_vertex, 1, MPI_UINT32_T, p, MPI_COMM_WORLD );
                to_weight_t received_to_weight;
                MPI_Status status;
                vertex_id_t cheapest_received_from = G->rank;
                for(int p = 0; p < G->nproc - 1; p++) {
                  MPI_Recv(&received_to_weight, 1,mpi_weight, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
                  if(received_to_weight.to != UINT32_MAX && received_to_weight.weight < cur_full_edge->weight){
                    cur_full_edge->to = received_to_weight.to;
                    cur_full_edge->weight = received_to_weight.weight;
                    cheapest_received_from = status.MPI_SOURCE;
                  }
                }
                vertex_id_t source = merge_ask(full_edges, cur_vertex, G, p, edge_parents, edges, changed, cheapest_received_from);
                if(source == G->rank){
                  added_edges[full_edges[i].edge] = true;
                }
              }
            }
            MPI_Bcast( &G->n, 1, MPI_UINT32_T, p, MPI_COMM_WORLD );
          }
          else{
            vertex_id_t component = 0;
            while(true){
              MPI_Bcast( &component, 1, MPI_UINT32_T, p, MPI_COMM_WORLD );
              if(component == G->n) break;
              to_weight_t to_weight;
              vertex_id_t index;
              vertex_id_t to;
              auto found = edges.find(component);
              if(found == edges.end() || found->second.local == UINT32_MAX){
                to_weight.to = UINT32_MAX;
              }else{
                index = found->second.local;
                to = full_edges[index].to;
                auto found_to = edges.find(to);
                if(found_to != edges.end() && found_to->second.edge != UINT32_MAX){
                  to = G->endV[find_in_edges(edge_parents,found_to->second.edge,G)];
                }
                to_weight.to = to;
                to_weight.weight = full_edges[index].weight;
              }
              MPI_Send(&to_weight, 1,mpi_weight,p, 0,MPI_COMM_WORLD);
              vertex_id_t source = merge_ask(full_edges, component, G, p, edge_parents, edges, changed);
              if(source == G->rank){
                added_edges[full_edges[index].edge] = true;
              }
            };

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
        auto found = components_map.find(foreign_index);
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
      return &trees;
    }
  }
}


extern "C" void* MST(graph_t *G) {
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