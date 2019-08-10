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
    fflush(stdout);

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
      printf("max_number = %u\n", max_number);
      printf("tree_mst.count = %u\n", trees_mst.size());

      vertex_id_t components_count = 0;
      for(vertex_id_t i = 0; i < trees_count; i++)
      {
        vertex_id_t component =  trees_mst[i][0];
        // printf("count=%d\n", trees_mst[i].size() - 1);
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
            // printf("count=%d first=%u\n", count, found->second);
            for(int k = 1; k < count; k++)
            {
              trees_mst[found->second].push_back(buf_edges[k]);
            }

          }


          // printf("count=%u from proc=%d\n", count, proc);
          // for(edge_id_t k= 0; k < count; k++){
          //   printf("[%u]=%u\n", k, buf_edges[k]);
          // }
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
              // printf("e[%u, %u]=%u\n", i, j - 1, trees_mst[i][j]);
              k++;
          }
      }

      trees_output->numTrees = trees_mst.size();
      trees_output->numEdges = number_of_edges;
      printf("numedges=%u\n",number_of_edges );
      printf("numTrees=%u\n",trees_output->numTrees );

    }else{
      MPI_Send(&trees_count, 1,MPI_UINT32_T,0, MPI_TAG, MPI_COMM_WORLD);
      for(vertex_id_t i = 0; i < trees_count; i++)
      {
          // vertex_id_t component = it->first;
          // vertex_id_t i = it->second;
          MPI_Send(&trees_mst[i][0], trees_mst[i].size(),MPI_UINT64_T, 0, MPI_TAG, MPI_COMM_WORLD);

      }
    }
    fflush(stdout);

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
      // if((from==60) ){
        // printf("from = %u gl_from_local=%u gl_to_local=%u to=%u %d\n", from,from_local+ first_vertex, to_local+ first_vertex, to, G->rank);
        // printf("!-!-!-!-!from = %u to=%u  %d\n", from,to, G->rank);
      //
      // }
      edge_id_t edge = find_in_edges(edge_parents, found->second.edge, G);
      if(to != G->endV[edge] ){
        // printf("from = %u gl_from_local=%u gl_to_local=%u to=%u %d\n", from,from_local+ first_vertex, to_local+ first_vertex, to, G->rank);
        // printf("!-!-!-!-! prev to = %u to=%u  %d\n", to,G->endV[edge], G->rank);

      }

      to = G->endV[edge];

      full_edges[local_from].to = to;

    }
    else{
      // printf("NOT FOUND %u rank = %d\n", from, G->rank);
      // fflush(stdout);
    }
    if(full_edges[local_from].weight == DBL_MAX){
      to = UINT32_MAX;
    }


    to_with_source[0] = to;
    to_with_source[1] = cheapest_received_from;

    MPI_Bcast( to_with_source, 2, MPI_UINT32_T, p, MPI_COMM_WORLD );
    if(from == 60 && from!=to){
      // printf("\nto = %u\n",to);
      // for(auto it = edges.cbegin(); it != edges.cend(); ++it)
      // {
      //     std::cout << it->first << " local=" << it->second.local << " edge=" << it->second.edge << "\n";
      // }
    }

  }
  else{
    MPI_Bcast( to_with_source, 2, MPI_UINT32_T, p, MPI_COMM_WORLD );
    to = to_with_source[0];
  }
  if(from == to || to == UINT32_MAX){
    // printf("DO not MERGE from=%u to=%u %d\n", from,to, G->rank);
    return UINT32_MAX;
  }
  changed = true;
  // printf("changed  %d\n", G->rank);
  // fflush(stdout);

  // if(G->rank == 0){
  //   printf("%u -> %u \n", from, to);
  // }

  auto found_from = edges.find(from);

  if(found_from == edges.end()){
    // printf("ERROR from=%u rank=%d\n",from,G->rank);
    // for(auto it = edges.cbegin(); it != edges.cend(); ++it)
    // {
    //     std::cout << it->first << " local=" << it->second.local << " edge=" << it->second.edge << "\n";
    // }
    // printf("RETURN rank = %d\n",G->rank);
    return UINT32_MAX;
  }
  // if((from == 44 || to == 44) && G->rank == 4){
  //   printf("!=!=!=from=%u to=%u\n", from,to);
  // }

  // if(p==G->rank){
  //   fflush(stdout);
  //   printf("%u -> %u, rank= %d\n", from,to, G->rank);
  //   fflush(stdout);
  //   if(from == 56)
  //     for(auto it = edges.cbegin(); it != edges.cend(); ++it)
  //     {
  //         std::cout << it->first << " local=" << it->second.local << " edge=" << it->second.edge << "\n";
  //     }
  //   fflush(stdout);
  // }

  // if(G->rank == 4 && from == 56){
  //   // fflush(stdout);
  //   // printf("%u -> %u, rank= %d\n", from,to, G->rank);
  //   // fflush(stdout);
  //   // if(from == 56)
  //     for(auto it = edges.cbegin(); it != edges.cend(); ++it)
  //     {
  //         std::cout << it->first << " local=" << it->second.local << " edge=" << it->second.edge << "\n";
  //     }
  //   fflush(stdout);
  // }




  fflush(stdout);
  // if((G->rank == 7 || G->rank == 5 ) && (from == 13 || from == 15)){
  //   printf("\nfrom = %u rank=%d\n",from,G->rank);
  //   for(auto it = edges.cbegin(); it != edges.cend(); ++it)
  //   {
  //       std::cout << it->first << " local=" << it->second.local << " edge=" << it->second.edge << "\n";
  //   }
  //
  //   // for(vertex_id_t j=0; j<G->local_m;j++){
  //   //
  //   // }
  //
  // }





  auto found_to = edges.find(to);

  if(found_to == edges.end()){
      // printf("NOT FOUND_TO %u  to %u %d received_from=%u\n", from,to, G->rank, to_with_source[1]);
    // if((from == 44 || to == 44) && G->rank == 4){
    //   printf("AAAAAAfrom=%u to=%u\n", from,to);
    // }

    loc_edge_t loc_edge;
    loc_edge.local = found_from->second.local;
    loc_edge.edge = found_from->second.edge;

    edges.insert( pair<vertex_id_t, loc_edge_t>(to, loc_edge) );
    if(loc_edge.edge != UINT64_MAX)
      G->endV[loc_edge.edge] = to;
    if(loc_edge.local != UINT32_MAX)
      full_edges[loc_edge.local].foreign = to;

  }else{
    // printf("FOUND_TO %u  to %u %d received_from=%u\n", from,to, G->rank, to_with_source[1]);

    // if((from == 44 || to == 44) && G->rank == 4){
    //   printf("AAAAAAfrom=%u to=%u\n", from,to);
    // }

    vertex_id_t from_local = found_from->second.local;
    vertex_id_t to_local = found_to->second.local;
    edge_id_t from_edge = found_from->second.edge;
    edge_id_t to_edge = found_to->second.edge;
    // if(false || G->rank==3){
      // for(auto it = edges.cbegin(); it != edges.cend(); ++it)
      // {
      //     std::cout << it->first << " local=" << it->second.local << " edge=" << it->second.edge << "\n";
      // }
      vertex_id_t printed_from = from_local == UINT32_MAX ? UINT32_MAX : from_local + first_vertex;
      vertex_id_t printed_to = to_local == UINT32_MAX ? UINT32_MAX : to_local + first_vertex;

    if(G->rank==4){
    // if((to == 60 || from==60 || (from_local + first_vertex == 24) || (from_local + first_vertex == 26) ) && G->rank==3){
      vertex_id_t printed_from = from_local == UINT32_MAX ? UINT32_MAX : from_local + first_vertex;
      vertex_id_t printed_to = to_local == UINT32_MAX ? UINT32_MAX : to_local + first_vertex;

      // printf("from = %u gl_from_local=%u gl_to_local=%u to=%u %d\n", from,printed_from, printed_to, to, G->rank);
      // exit(0);
    }

    // if(from_local == UINT32_MAX){
    //   if(to_local != UINT32_MAX){
    //     full_edges[to_local].parent = to_local;
    //
    //   }
    // }

    if(from_local != UINT32_MAX){
      if(to_local == UINT32_MAX){
        // printf("to_local nil %d\n", G->rank);
        full_edges[from_local].foreign = to;
        found_to->second.local = from_local;

      }else{
        // printf("to_local present from_local=%u to_local=%u %d\n", from_local, to_local, G->rank);
        full_edges[from_local].parent = to_local;
        // if(G->rank==3){
        //   printf("from = %u gl_from_local=%u gl_to_local=%u to=%u %d\n", from,from_local+ first_vertex, to_local+ first_vertex, to, G->rank);
        // }

        full_edges[from_local].foreign = UINT32_MAX;
        if(full_edges[from_local].weight < full_edges[to_local].weight){
          full_edges[to_local].to = full_edges[from_local].to;
          full_edges[to_local].weight = full_edges[from_local].weight;
          full_edges[to_local].edge = full_edges[from_local].edge;
          // printf("!!!vertex %u of  %u |  %u edge=%u\n", VERTEX_TO_GLOBAL(i, G->n, G->nproc, G->rank),
          //  global_first_component_index,second_component_index, j);
        }
      }
    }

    if(from_edge == UINT64_MAX){
      if(to_edge == UINT64_MAX){
      }else{
        found_from->second.edge = to_edge;
        // found_to->second.edge = to_edge;

        // edge_parents[from_edge] = found_to->second.edge;
      }
    }else{
      if(to_edge == UINT64_MAX){
        //local есть
        G->endV[from_edge] = to;
        found_to->second.edge = from_edge;

      }else{
        // printf("!!!!\n");
        G->endV[from_edge] = to;
        found_to->second.edge = to_edge;
        edge_parents[from_edge] = found_to->second.edge;

      }
    }

    // printf("RANK=%d  found_from_first=%u found_to_first=%u found_from=%u found_to=%u\n",
    // G->rank, found_from->first, found_to->first, found_from->second, found_to->second);
    // edge_parents[found_from->second] = found_to->second;
  }
  // printf("END=%d\n", G->rank);
  // fflush(stdout);
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
vertex_id_t iter = 0;

extern "C" void* MST_boruvka(graph_t *G) {
    int rank = G->rank, size = G->nproc;
    bool changed;
    full_edge_t *full_edges = new full_edge_t[G->local_n];
    edge_id_t *edge_parents = new edge_id_t[G->local_m];
    map <vertex_id_t, loc_edge_t> edges;
    bool *added_edges = new bool[G->local_m];
    memset(added_edges,0,G->local_m * sizeof(bool));

    // map <vertex_id_t,vertex_id_t> components_map;
    // print_source_graph(G);
    // exit(0);


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
    // for(vertex_id_t i = 0; i < G->local_n; i++) {
    //   for (edge_id_t j = G->rowsIndices[i]; j < G->rowsIndices[i+1]; j++) {
    //     if(first_vertex + i == 0){
    //       G->endV[j] = 2;
    //     }
    //
    //     if(first_vertex + i == 2){
    //       G->endV[j] = 0;
    //     }
    //
    //     if(first_vertex + i == 1){
    //       G->endV[j] = 3;
    //     }
    //     if(first_vertex + i == 3){
    //       G->endV[j] = 1;
    //     }
    //
    //
    //     // if((first_vertex + i == 1 || first_vertex + i == 3) &&  (G->endV[j] == 0 || G->endV[j] == 2)){
    //     //   G->endV[j] = first_vertex + i;
    //     // }
    //
    //   }
    // }

    // print_source_graph(G);
    //
    // MPI_Barrier(MPI_COMM_WORLD);
    // fflush(stdout);
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
    printf("\n");
    // if(G->rank == 7)
    //   for(auto it = edges.cbegin(); it != edges.cend(); ++it)
    //   {
    //       std::cout << it->first << " local=" << it->second.local << " edge=" << it->second.edge << "\n";
    //   }




    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    while(true){
        iter++;
        // if(G->rank == 0)
          // printf("iter %u\n", iter);
        fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);
        fflush(stdout);
        // if(iter == 10) {
        //   fflush(stdout);
        //   exit(0);
        // }
        bool changed = false;
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
            // vertex_id_t global_first_component_index = component_from_find(full_edges, first_component_index, first_vertex, G, edge_parents, edges);
            // vertex_id_t global_first_component_index = full_edges[i].foreign == UINT32_MAX ? first_component_index + first_vertex : full_edges[i].foreign;
            vertex_id_t global_first_component_index = full_edges[first_component_index].foreign;

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
              // printf("!!!vertex %u of  %u |  %u edge=%u\n", VERTEX_TO_GLOBAL(i, G->n, G->nproc, G->rank),
              //  global_first_component_index,second_component_index, j);
            }
            // assign_cheapest(full_edges, G, i, j);
            // foreign
            // ребро своё
            //
            //
          }
          // printf("full_edges[%u, p = %u] = %u\n", VERTEX_TO_GLOBAL(i, G->n, G->nproc, G->rank),
          // component_from_find(full_edges, i, first_vertex, G, edge_parents, edges), full_edges[i].to);
        }
        fflush(stdout);

        for(int p = 0; p < G->nproc; p++) {

          if(p == G->rank){
            for(vertex_id_t i = 0; i < G->local_n; i++) {
              // if(iter==3){
                // printf("RANK= %d trying %u\n", G->rank, first_vertex + i);
                // fflush(stdout);
              // }
              full_edge_t *cur_full_edge = &full_edges[i];
              if(cur_full_edge->foreign !=UINT32_MAX && cur_full_edge->foreign == i + first_vertex){
                // printf("-------trying %u RANK= %d\n", first_vertex + i, G->rank);
                // fflush(stdout);

                vertex_id_t cur_vertex = first_vertex + i;
                // printf("I WANT %u foreign=%u parent=%u\n", cur_vertex, cur_full_edge->foreign, first_vertex + cur_full_edge->parent);
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
                // printf("%u results_to=%u received_from=%u my_rank=%d\n", i+ first_vertex, cur_full_edge->to, cheapest_received_from, G->rank);

                // printf("rank=%d my=%u merge with=%u\n", G->rank,VERTEX_TO_GLOBAL(i, G->n, G->nproc, G->rank),cur_full_edge->to);
                vertex_id_t source = merge_ask(full_edges, cur_vertex, G, p, edge_parents, edges, changed, cheapest_received_from);
                // printf("source=%u rank=%d OWNER\n", source, G->rank);
                if(source == G->rank){
                  added_edges[full_edges[i].edge] = true;
                  // printf("%")
                  fflush(stdout);
                  // printf("MERGE %u TO %u | edge=%u weight=%lf rank=%d\n",
                  // i + first_vertex, full_edges[i].to,
                  // edge_to_global(full_edges[i].edge,G), full_edges[i].weight, G->rank);
                  fflush(stdout);
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
              if(component == G->n) break;
              to_weight_t to_weight;
              vertex_id_t index;
              vertex_id_t to;
              auto found = edges.find(component);
              if(found == edges.end() || found->second.local == UINT32_MAX){
                // printf("NOT FOUND=%u %d\n", component, G->rank );
                // fflush(stdout);

                to_weight.to = UINT32_MAX;
              }else{
                // printf("FOUND=%u %d\n", component, G->rank );
                // fflush(stdout);
                index = found->second.local;

                // vertex_id_t to = first_vertex +  find(full_edges, index, G);
                // if(to == from)
                //   to_weight.to = UINT32_MAX;
                // else
                to = full_edges[index].to;
                auto found_to = edges.find(to);
                if(found_to != edges.end() && found_to->second.edge != UINT32_MAX){
                  to = G->endV[find_in_edges(edge_parents,found_to->second.edge,G)];
                }
                if(full_edges[index].to != to){
                  // printf("CHANGED FROM=%u TO=%u\n",full_edges[index].to,to );
                }
                to_weight.to = to;
                // full_edges[index].to = to;

                // to_weight.to =  full_edges[index].to;
                if(full_edges[index].parent != index){
                  printf("ERROR!!!!!!!!");
                  exit(0);
                }
                to_weight.weight = full_edges[index].weight;
              }
              MPI_Send(&to_weight, 1,mpi_weight,p, 0,MPI_COMM_WORLD);
              vertex_id_t source = merge_ask(full_edges, component, G, p, edge_parents, edges, changed);
              // printf("source=%u rank=%d\n", source, G->rank);
              if(source == G->rank){
                added_edges[full_edges[index].edge] = true;
                // full_edges[index].added_edge = full_edges[index].edge;

                // if(component == 60){
                //   printf("\component = %u rank=%d\n",component, G->rank);
                //   printf("\n");
                //   // for(vertex_id_t k = 0; k < G->local_n; k++){
                //   //   printf("[%u].parent=%u\n",k+first_vertex, full_edges[k].parent + first_vertex);
                //   // }
                //   for(auto it = edges.cbegin(); it != edges.cend(); ++it)
                //   {
                //       std::cout << it->first << " local=" << it->second.local << " edge=" << it->second.edge << "\n";
                //   }
                //   fflush(stdout);
                // }
                // fflush(stdout);
                // if(component == 56){
                //   printf("\component = %u rank=%d\n",component, G->rank);
                //   printf("\n");
                //   // for(vertex_id_t k = 0; k < G->local_n; k++){
                //   //   printf("[%u].parent=%u\n",k+first_vertex, full_edges[k].parent + first_vertex);
                //   // }
                //   for(auto it = edges.cbegin(); it != edges.cend(); ++it)
                //   {
                //       std::cout << it->first << " local=" << it->second.local << " edge=" << it->second.edge << "\n";
                //   }
                //   fflush(stdout);
                // }
                // printf("--------MERGE_CHILD %u TO %u | edge=%u of vertex %u foreign=%u  par=%u weight=%lf rank=%d \n", component,
                // to, edge_to_global(full_edges[index].edge,G), first_vertex + index,
                // full_edges[index].foreign,first_vertex + full_edges[index].parent, full_edges[index].weight, G->rank);
                fflush(stdout);

              }
            };

          }
          MPI_Barrier(MPI_COMM_WORLD);
        }
        // exit(0);

        MPI_Barrier(MPI_COMM_WORLD);

    if(!changed){
      // printf("CHANGED \n ");
      // fflush(stdout);
      // exit(0);

      // MPI_Barrier(MPI_COMM_WORLD);
      // exit(0);
      map <vertex_id_t,vertex_id_t> components_map;
      trees.clear();
      vertex_id_t count = 0;
      for(vertex_id_t i = 0; i < G->local_n;i++){
        vertex_id_t global_vertex = i + first_vertex;
        vertex_id_t component_index =  find(full_edges, i, G);
        // edge_id_t edges_count = G->rowsIndices[component_index+1] - G->rowsIndices[component_index];
        vertex_id_t foreign_index =  full_edges[component_index].foreign;
        // if(edges_count == 0){
        //   printf("AAAAA\n");
        //   components_map.insert(pair<vertex_id_t, vertex_id_t>(foreign_index, count));
        //   trees.push_back(vector<edge_id_t>());
        //   trees[count].push_back(foreign_index);
        //   count++;
        //   continue;
        // }

        auto found = components_map.find(foreign_index);
        vertex_id_t address = 0;

        if(found == components_map.end()){
          // printf("ADDED = %u rank=%d\n", foreign_index, rank);
          components_map.insert(pair<vertex_id_t, vertex_id_t>(foreign_index, count));
          trees.push_back(vector<edge_id_t>());
          trees[count].push_back(foreign_index);
          address = count;
          // trees[count].push_back(edge_to_global(edge,G));
          count++;
        }else{
          address = found->second;

          // trees[found->second].push_back(edge_to_global(edge,G));
        }



        for(edge_id_t j = G->rowsIndices[i]; j < G->rowsIndices[i+1];j++){
          if(!added_edges[j]){
            continue;
          }
          trees[address].push_back(edge_to_global(j,G));
        }
      }

      // printf("TREES_SIZE=%u rank=%d\n", trees.size(), G->rank);

      for(vertex_id_t i = 0; i < trees.size();i++){
        for(vertex_id_t j = 0; j < trees[i].size();j++){
          // printf("[%u,%u]=%u rank=%d\n", i, j, trees[i][j], G->rank);
        }
      }

      // printf("END rank=%d \n", G->rank);



      // for(vertex_id_t i = 0; i < G->local_n;i++){
      //   if(i + first_vertex != component) continue;
      //   trees.push_back(vector<edge_id_t>());
      //   components_map.insert(pair <vertex_id_t, vertex_id_t>(component,trees.size() - 1));
      // }
      // printf("CHANGED SECOND \n ");
      // MPI_Barrier(MPI_COMM_WORLD);
      //
      // trees.push_back(vector<edge_id_t>());

      // for(edge_id_t i = 0; i < G->local_m; i++){
      //   edge_id_t edge = i;
      //   if(!added_edges[edge]){
      //     continue;
      //   }
      //   // printf("CHANGED edge=%u rank=%d\n ", i, G->rank);
      //   // vertex_id_t component = find(full_edges, G->endV[edge], G);
      //   // auto addr = components_map.find(component);
      //   trees[0].push_back(edge_to_global(edge,G));
      //   // trees[addr->second].push_back(edge_to_global(edge,G));
      // }
      // printf("END rank=%d\n ",  G->rank);

      free(full_edges);
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
