#include <vector>
#include <mpi.h>
#include <float.h>
#include <stdlib.h>
#include <assert.h>
#include "defs.h"
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <omp.h>
#include <unistd.h>

#include <map>

#define UINT32_MAX  (0xffffffff)
#define UINT64_MAX  (0xffffffffffffffff)


using namespace std;
map <vertex_id_t,vertex_id_t> foreign_vertexes;


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
            printf("%d to %d, weight= %lf, edge_number= %d  \n", a, G->endV[j], G-> weights[j], j);
        }

      }
      // MPI_Barrier(MPI_COMM_WORLD);
      printf("\n");
      // }
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
    result_t &trees_mst = *reinterpret_cast<result_t*>(result);
    int size, rank;
    int tree_size = 0;
    int MPI_TAG = 99;
    edge_id_t edge_buf;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    for(int tree_count = 0; (unsigned int)tree_count < trees_mst.size(); ++tree_count) //trees loop
        for(int proc_count = 1; proc_count < size; ++proc_count)         //processes loop
            if(rank == proc_count || rank == 0)
            {
                MPI_TAG = tree_count;
                //send to 0-process number of edges in tree
                if (rank == proc_count)
                {
                    tree_size = trees_mst[tree_count].size();
                    MPI_Send(&tree_size, 1, MPI_INT, 0, MPI_TAG, MPI_COMM_WORLD);
                }
                // 0-process recieve tree_size from any process
                else
                {
                    MPI_Recv(&tree_size, 1, MPI_INT, proc_count, MPI_TAG, MPI_COMM_WORLD, &status);
                }
                // 0-process recieve edge_id from any process
                for(int edge_count = 0; edge_count < tree_size; ++edge_count) //edges loop
                {
                    MPI_TAG = edge_count;
                    if(rank == proc_count)
                    {
                         MPI_Send(&(trees_mst[tree_count])[edge_count], 1, MPI_UINT64_T, 0, MPI_TAG, MPI_COMM_WORLD);
                    }
                    else
                    {
                        MPI_Recv(&edge_buf, 1, MPI_UINT64_T, proc_count, MPI_TAG, MPI_COMM_WORLD, &status);
                        trees_mst[tree_count].push_back(edge_buf);
                    }
                }
            }
    if (rank == 0)
    {
        trees_output->p_edge_list = (edge_id_t *)malloc(trees_mst.size()*2 * sizeof(edge_id_t));
        edge_id_t number_of_edges = 0;
        for (vertex_id_t i = 0; i < trees_mst.size(); i++) number_of_edges += trees_mst[i].size();
        trees_output->edge_id = (edge_id_t *)malloc(number_of_edges * sizeof(edge_id_t));
        trees_output->p_edge_list[0] = 0;
        trees_output->p_edge_list[1] = trees_mst[0].size();
        for (vertex_id_t i = 1; i < trees_mst.size(); i++) {
            trees_output->p_edge_list[2*i] = trees_output->p_edge_list[2*i-1];
            trees_output->p_edge_list[2*i +1] = trees_output->p_edge_list[2*i-1] + trees_mst[i].size();
        }
        int k = 0;
        for (vertex_id_t i = 0; i < trees_mst.size(); i++) {
            for (edge_id_t j = 0; j < trees_mst[i].size(); j++) {
                trees_output->edge_id[k] = trees_mst[i][j];
                // printf("e[%u, %u]=%u\n", i, j, trees_mst[i][j]);
                k++;
            }
        }

        trees_output->numTrees = trees_mst.size();
        trees_output->numEdges = number_of_edges;
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

// vertex_id_t find(full_edge_t* full_edges, vertex_id_t i, graph_t* G)
// {
//   vertex_id_t new_i;
//   // printf("full_edges[%u].end=%u\n", i, full_edges[i].end);
//
//   while(!full_edges[i].end){
//     // printf("i= %u\n", i);
//     new_i = full_edges[i].parent;
//     full_edges[i].parent = full_edges[new_i].parent;
//     i = full_edges[new_i].parent;
//   }
//   // printf("full_edges[%u].parent=%u\n", i, full_edges[i].parent);
//   // fflush(stdout);
//   return i;
//
//   // return full_edges[i].parent;
//   // return VERTEX_OWNER(i, G->n, G->nproc);
// }

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
               edge_id_t *edge_parents, map <vertex_id_t, edge_id_t> &edges,
               bool &changed, vertex_id_t cheapest_received_from = 0){
  // vertex_id_t from = i;
  // vertex_id_t from = i;
  vertex_id_t first_vertex = VERTEX_TO_GLOBAL(0, G->n, G->nproc, G->rank);
  vertex_id_t to;
  vertex_id_t to_with_source[2];
  // MPI_Barrier(MPI_COMM_WORLD);
  // exit(0);
  if(p == G->rank){
    vertex_id_t local_from = from - first_vertex;
    to = full_edges[local_from].to;
    auto found = edges.find(to);

    // edge_id_t edge = find_in_edges(edge_parents, found->second, G);
    // to = G->endV[edge];
    if(found != edges.end()){
      edge_id_t edge = find_in_edges(edge_parents, found->second, G);
      to = G->endV[edge];
    }
    else{
      printf("NOT FOUND %u rank = %d\n", from, G->rank);
      fflush(stdout);
    }

    to_with_source[0] = to;
    to_with_source[1] = cheapest_received_from;

    MPI_Bcast( to_with_source, 2, MPI_UINT32_T, p, MPI_COMM_WORLD );
    if(VERTEX_OWNER(to, G->n, G->nproc) == p){
      // edge_id_t edge = find_in_edges(edge_parents, edges.find(to)->second, G);
      // to = G->endV[edge];

      full_edges[local_from].parent = to;
    }
    else{
      full_edges[local_from].foreign = to;
    }
  }
  else{
    MPI_Bcast( to_with_source, 2, MPI_UINT32_T, p, MPI_COMM_WORLD );
    to = to_with_source[0];
  }
  if(from == to || to == UINT32_MAX){
    return UINT32_MAX;
  }
  changed = true;
  if(G->rank == 0){
    printf("MERGE %u TO %u | received_from=%u \n", from, to, to_with_source[1]);
  }

  auto found_from = edges.find(from);
  auto found_to = edges.find(to);
  if(found_from == edges.end()){
    return UINT32_MAX;
  }

  if(found_to == edges.end()){
    edges.insert( pair<vertex_id_t,edge_id_t>(found_to->first, found_from->second) );
    printf("DDDDDDDDDDD=%d\n", G->rank);

    G->endV[found_from->second] = found_to->first;
  }else{
    // printf("RANK=%d  found_from_first=%u found_to_first=%u found_from=%u found_to=%u\n",
    // G->rank, found_from->first, found_to->first, found_from->second, found_to->second);
    edge_parents[found_from->second] = found_to->second;
  }
  return to_with_source[1];
  // edges.erase(found_from);
}


// void merge_answer(int p){
//   MPI_Bcast( &to, 1, MPI_UINT32_T, p, MPI_COMM_WORLD );
//
// }


// void forward_ffff(){
//
// }

extern "C" void* MST_boruvka(graph_t *G) {
    int rank = G->rank, size = G->nproc;
    bool changed;
    full_edge_t *full_edges = new full_edge_t[G->local_n];
    edge_id_t *edge_parents = new edge_id_t[G->local_m];
    bool *added_edges = new bool[G->local_m];
    memset(added_edges,0,G->local_m * sizeof(bool));
    map <vertex_id_t, edge_id_t> edges;

    // print_source_graph(G);

    int weight_nitems = 2;
    int weight_blocklengths[weight_nitems] = {1,1};
    MPI_Datatype weight_types[weight_nitems] = {MPI_UINT32_T, MPI_DOUBLE};
    MPI_Datatype mpi_weight;
    MPI_Aint     weight_offsets[weight_nitems];

    weight_offsets[0] = offsetof(to_weight_t, to);
    // weight_offsets[1] = offsetof(to_weight_t, rank);
    weight_offsets[1] = offsetof(to_weight_t, weight);

    MPI_Type_create_struct(weight_nitems, weight_blocklengths, weight_offsets, weight_types, &mpi_weight);
    MPI_Type_commit(&mpi_weight);



    vertex_id_t first_vertex = VERTEX_TO_GLOBAL(0, G->n, G->nproc, G->rank);

    for(vertex_id_t i = 0; i < G->local_n; i++){
      full_edges[i].rank = 0;
      full_edges[i].parent = i;
      full_edges[i].foreign = UINT32_MAX;
    }

    for(edge_id_t i = 0; i < G->local_m; i++){
      vertex_id_t vertex = G->endV[i];
      // if(VERTEX_OWNER(vertex, G->n, G->nproc) == rank){
      //   // foreign_edges[i].parent = i;
      //   foreign_edges[i].parent = i;
      //   foreign_edges[i].end = true;
      //   continue;
      // }
      auto found = edges.find(vertex);
      if(found == edges.end()){
        edges.insert( pair<vertex_id_t,edge_id_t>(vertex,i) );
        edge_parents[i] = i;
        // printf("new %i", )

      }else{
        edge_parents[i] = found->second;
      }
    }


    vertex_id_t iter = 0;
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    while(true){
        if(G->rank == 0)
          printf("iter %u\n", iter);
        fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);
        fflush(stdout);
        if(iter == 10) {
          fflush(stdout);
          exit(0);
        }
        bool changed = false;
        iter++;
        for(vertex_id_t i = 0; i < G->local_n; i++){
          full_edges[i].to = UINT32_MAX;
          full_edges[i].weight = DBL_MAX;

        }
        // printf("rank=%d REACHED\n", rank);
        // fflush(stdout);
        // MPI_Barrier(MPI_COMM_WORLD);
        // if(iter == 2) exit(1);

        // memset(count_to,0,size * sizeof(vertex_id_t));
        for(vertex_id_t i = 0; i < G->local_n; i++) {
          for (edge_id_t j = G->rowsIndices[i]; j < G->rowsIndices[i+1]; j++) {

            vertex_id_t first_component_index =  find(full_edges, i, G);
            vertex_id_t global_first_component_index = component_from_find(full_edges, first_component_index, first_vertex, G, edge_parents, edges);
            // vertex_id_t second_vertex = G->endV[j];
            vertex_id_t second_component_index = G->endV[find_in_edges(edge_parents,j,G)];

            // printf("BEGIN = %u |  %u %u\n", first_component_index, second_component_index, j);
            // printf("rank=%d REACHED\n", rank);
            // fflush(stdout);
            if(first_component_index == UINT32_MAX){
              // printf("MAX\n");
            }

            if(global_first_component_index == second_component_index){
              // printf("%u = %u FROM RANK %d\n",  VERTEX_TO_GLOBAL(i, G->n, G->nproc, G->rank), second_component_index, G->rank);
              // fflush(stdout);
              continue;
            }

            // if ( (foreign == UINT32_MAX ? first_vertex + first_component_index : foreign)  == second_component_index){
            //   continue;
            // }
            if(full_edges[first_component_index].weight > G->weights[j]){
              full_edges[first_component_index].to = second_component_index;
              full_edges[first_component_index].weight = G->weights[j];
              full_edges[first_component_index].edge = j;
              // printf("!!!first_component_index = %u |  %u %u\n", VERTEX_TO_GLOBAL(i, G->n, G->nproc, G->rank), second_component_index, j);
            }
            // assign_cheapest(full_edges, G, i, j);
            // foreign
            // ребро своё
            //
            //
          }
          printf("full_edges[%u, p = %u] = %u\n", VERTEX_TO_GLOBAL(i, G->n, G->nproc, G->rank),
          component_from_find(full_edges, i, first_vertex, G, edge_parents, edges), full_edges[i].to);
        }
        fflush(stdout);


        for(int p = 0; p < G->nproc; p++) {
          if(p == G->rank){
            for(vertex_id_t i = 0; i < G->local_n; i++) {
              // printf("RANK= %d trying %u\n", G->rank, first_vertex + i);
              full_edge_t *cur_full_edge = &full_edges[i];
              if(cur_full_edge->foreign == UINT32_MAX){
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
                    // cur_full_edge->foreign_rank = received_to_weight.rank;
                  }
                }
                // printf("rank=%d my=%u merge with=%u\n", G->rank,VERTEX_TO_GLOBAL(i, G->n, G->nproc, G->rank),cur_full_edge->to);
                vertex_id_t source = merge_ask(full_edges, cur_vertex, G, p, edge_parents, edges, changed, cheapest_received_from);
                if(source == G->rank){
                  added_edges[full_edges[i].edge] = true;
                }

                //
                // MPI_Bcast( &cur_full_edge->parent, 1, MPI_UINT32_T, p, MPI_COMM_WORLD );
              }
            }
            MPI_Bcast( &G->n, 1, MPI_UINT32_T, p, MPI_COMM_WORLD );
          }
          else{
            vertex_id_t component = 0;
            // printf("HERE=%d\n", G->rank );
            // fflush(stdout);

            while(true){
              MPI_Bcast( &component, 1, MPI_UINT32_T, p, MPI_COMM_WORLD );
              // printf("HERE=%d\n", G->rank );
              // fflush(stdout);
              if(component == G->n) break;
              to_weight_t to_weight;
              vertex_id_t index = UINT32_MAX;
              for(vertex_id_t i = 0; i < G->local_n; i++) {
                // if(full_edges[i].foreign == component)
                // printf("comp1 %d\n", G->rank);
                if(component_from_find(full_edges, i, first_vertex, G, edge_parents, edges) == component)
                {
                  index = i;
                }
                // printf("comp2 %d\n", G->rank);

              }
              edge_id_t edge;
              if(index == UINT32_MAX){
                // printf("NOT FOUND rank=%d NOT FOUND for component=%u\n", G->rank, component );
                // fflush(stdout);

                to_weight.to = UINT32_MAX;
              }else{
                // vertex_id_t home = find_in_edges(foreign_edges, found->second, G).home_reference;
                // full_edges[home]
                // vertex_id_t index = found -> first;
                // printf("FOUND rank=%d to=%u for component=%u\n", G->rank,full_edges[index].to, component);
                fflush(stdout);

                vertex_id_t to = full_edges[index].to;
                edge = find_in_edges(edge_parents, edges.find(to)->second, G);
                to = G->endV[edge];
                if(to == component){
                  // printf("OLD rank=%d to=%u for component=%u\n", G->rank, to, component);
                  to_weight.to = UINT32_MAX;
                }else{
                  // printf("FOUND rank=%d to=%u for component=%u\n", G->rank, to, component);
                  to_weight.to =  full_edges[index].to;//first_vertex + home;//full_edges[home].to;
                  to_weight.weight = full_edges[index].weight;
                }
                // to_weight.rank = full_edges[index].rank;

              }
              MPI_Send(&to_weight, 1,mpi_weight,p, 0,MPI_COMM_WORLD);
              // merge_answer(component, p);

              vertex_id_t source = merge_ask(full_edges, component, G, p, edge_parents, edges, changed);
              if(source == G->rank){
                added_edges[edge] = true;
              }


              // for(vertex_id_t i = 0; i < G->local_n; i++) {
              //
              // }
            };

          }
          MPI_Barrier(MPI_COMM_WORLD);
        }
        // exit(0);

        MPI_Barrier(MPI_COMM_WORLD);


    map <vertex_id_t,vertex_id_t> components_map;


    if(!changed){
      trees.clear();
      printf("CHANGED \n ");
      MPI_Barrier(MPI_COMM_WORLD);

      // for(vertex_id_t i = 0; i < G->local_n;i++){
      //   if(i + first_vertex != component) continue;
      //   trees.push_back(vector<edge_id_t>());
      //   components_map.insert(pair <vertex_id_t, vertex_id_t>(component,trees.size() - 1));
      // }
      // printf("CHANGED SECOND \n ");
      // MPI_Barrier(MPI_COMM_WORLD);
      //
      trees.push_back(vector<edge_id_t>());

      for(edge_id_t i = 0; i < G->local_m; i++){
        edge_id_t edge = i;
        if(!added_edges[edge]){
          continue;
        }
        // printf("CHANGED edge=%u rank=%d\n ", i, G->rank);
        // vertex_id_t component = find(full_edges, G->endV[edge], G);
        // auto addr = components_map.find(component);
        trees[0].push_back(edge_to_global(edge,G));
        // trees[addr->second].push_back(edge_to_global(edge,G));
      }
      printf("END rank=%d\n ",  G->rank);

      free(full_edges);
      free(added_edges);
      fflush(stdout);
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
