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
#define MPI_UNITE_REQUEST  1
#define MPI_UNITE_RESPONSE 2


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

vertex_id_t find(full_edge_t* full_edges, vertex_id_t i, graph_t* G)
{
  vertex_id_t new_i;
  // printf("full_edges[%u].end=%u\n", i, full_edges[i].end);

  while(!full_edges[i].end){
    // printf("i= %u\n", i);
    new_i = full_edges[i].parent;
    full_edges[i].parent = full_edges[new_i].parent;
    i = full_edges[new_i].parent;
  }
  // printf("full_edges[%u].parent=%u\n", i, full_edges[i].parent);
  // fflush(stdout);
  return i;

  // return full_edges[i].parent;
  // return VERTEX_OWNER(i, G->n, G->nproc);
}

vertex_id_t find_in_edges(foreign_edges_t* foreign_edges, vertex_id_t i, graph_t* G)
{
  vertex_id_t new_i;
  while(i != foreign_edges[i].parent){
    new_i = foreign_edges[i].parent;
    foreign_edges[i].parent = foreign_edges[new_i].parent;
    i = foreign_edges[new_i].parent;
  }
  return G->endV[i];
}

extern "C" void* MST_boruvka(graph_t *G) {
    int rank = G->rank, size = G->nproc;
    bool changed;
    full_edge_t *full_edges = new full_edge_t[G->local_n];
    bool *added_edges = new bool[G->local_m];
    memset(added_edges,0,G->local_m * sizeof(bool));
    map <vertex_id_t, edge_id_t> foreign_components;
    foreign_edges_t *foreign_edges = new foreign_edges_t[G->local_m];

    from_to_t *all = new from_to_t[size];
    from_to_t *recv_all = new from_to_t[size];
    from_to_t *all_result = new from_to_t[size];
    vertex_id_t *send_to = new vertex_id_t[size];
    MPI_Request *requests = new MPI_Request[G->local_n];
    MPI_Status *statuses = new MPI_Status[G->local_n];

    vertex_id_t receive_count;
    int *recounts = new int[size];
    for(int i=0; i < size; i++ ){
      recounts[i] = 1;
    }

    print_source_graph(G);


    int nitems = 2;
    int blocklengths[nitems] = {1,1};
    MPI_Datatype types[nitems] = {MPI_UINT32_T, MPI_UINT32_T};
    MPI_Datatype mpi_edge;
    MPI_Aint     offsets[nitems];

    offsets[0] = offsetof(from_to_t, from);
    offsets[1] = offsetof(from_to_t, to);
    // offsets[2] = offsetof(rank_comparator_t, rank);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_edge);
    MPI_Type_commit(&mpi_edge);

    int weight_nitems = 3;
    int weight_blocklengths[weight_nitems] = {1,1,1};
    MPI_Datatype weight_types[weight_nitems] = {MPI_UINT32_T, MPI_UINT32_T, MPI_DOUBLE};
    MPI_Datatype mpi_weight;
    MPI_Aint     weight_offsets[weight_nitems];

    weight_offsets[0] = offsetof(from_to_with_weight_t, from);
    weight_offsets[1] = offsetof(from_to_with_weight_t, to);
    weight_offsets[2] = offsetof(from_to_with_weight_t, weight);

    MPI_Type_create_struct(weight_nitems, weight_blocklengths, weight_offsets, weight_types, &mpi_weight);
    MPI_Type_commit(&mpi_weight);



    vertex_id_t first_vertex = VERTEX_TO_GLOBAL(0, G->n, G->nproc, G->rank);

    for(vertex_id_t i = 0; i < G->local_n; i++){
      full_edges[i].parent =  first_vertex + i;
      full_edges[i].end =  true;
      // full_edges[i].rank = 0;
    }

    for(edge_id_t i = 0; i < G->local_m; i++){
      vertex_id_t vertex = G->endV[i];
      if(VERTEX_OWNER(vertex, G->n, G->nproc) == rank){
        // foreign_edges[i].parent = i;
        continue;
      }
      auto  found = foreign_components.find(vertex);
      if(found == foreign_components.end()){
        foreign_components.insert( pair<vertex_id_t,edge_id_t>(vertex,i) );
        foreign_edges[i].parent = i;
        foreign_edges[i].home_reference = UINT32_MAX;

      }else{
        foreign_edges[i].parent = found->second;
      }
    }


    vertex_id_t iter  = 0;

    while(true){
        bool changed = false;
        iter++;
        for(vertex_id_t i = 0; i < G->local_n; i++){
          // full_edges[i].to = UINT32_MAX;
          full_edges[i].weight = DBL_MAX;

        }
        // printf("rank=%d REACHED\n", rank);
        // fflush(stdout);
        // MPI_Barrier(MPI_COMM_WORLD);
        // if(iter == 2) exit(1);

        // memset(count_to,0,size * sizeof(vertex_id_t));
        for(int p = 0; p < size; p++){
        if(p == rank){
          // printf("rank = %d\n",rank );
        for(vertex_id_t i = 0; i < G->local_n; i++) {
          for (edge_id_t j = G->rowsIndices[i]; j < G->rowsIndices[i+1]; j++) {

            vertex_id_t local_first =  find(full_edges, i, G);
            vertex_id_t first_component_index =   full_edges[local_first].parent;
            vertex_id_t second_component_index;
            vertex_id_t vertex = G->endV[j];
            if(vertex >= first_vertex && vertex < first_vertex + G->local_n )
            second_component_index =  full_edges[find(full_edges, vertex - first_vertex, G)].parent;
            else
            second_component_index =   find_in_edges(foreign_edges, j, G); //G->endV[j]; //compute_second_component_index(full_edges, G, j);
            // printf("BEGIN = %u |  %u %u\n", first_component_index, second_component_index, j);
            // printf("rank=%d REACHED\n", rank);
            // fflush(stdout);


            if (first_component_index == second_component_index){
              continue;
            }
            // first_component_index = VERTEX_TO_GLOBAL(i, G->n, G->nproc, G->rank) - VERTEX_TO_GLOBAL(0, G->n, G->nproc, G->rank);
            // VERTEX_TO_GLOBAL(i, G->n, G->nproc, G->rank)
            // printf("first_component_index = %u |  %u %u\n", first_component_index, second_component_index, j);
            // edge_id_t edge = full_edges[local_first].edge;

            if (full_edges[local_first].weight > G->weights[j]){
              // printf("first_component_index = %u |  %u %u\n", VERTEX_TO_GLOBAL(i, G->n, G->nproc, G->rank), second_component_index, j);
              full_edges[local_first].to = second_component_index;
              full_edges[local_first].weight = G->weights[j];

              // printf("!!!first_component_index = %u |  %u %u\n", VERTEX_TO_GLOBAL(i, G->n, G->nproc, G->rank), second_component_index, j);
            }
            fflush(stdout);
            // assign_cheapest(full_edges, G, i, j);
          }
        }
      }
        MPI_Barrier(MPI_COMM_WORLD);
      }


      // MPI_Barrier(MPI_COMM_WORLD);
      memset(send_to,0,size * sizeof(vertex_id_t));
      // MPI_Barrier(MPI_COMM_WORLD);


      for(vertex_id_t i = 0; i < G->local_n; i++) {

        // // printf("local_form begin\n");
        // // fflush(stdout);
        // MPI_Barrier(MPI_COMM_WORLD);

        vertex_id_t local_from  =  find(full_edges,i,G);


        full_edge_t *edge = &full_edges[local_from];

        vertex_id_t p = edge->parent;
        vertex_id_t to = edge->to;
        weight_t weight = edge->weight;

        int owner = VERTEX_OWNER(p, G->n, G->nproc);
        if(owner == rank) continue;


        // printf("rank = %d owner = %d local_from = %d \n",rank, owner, local_from);
        // fflush(stdout);

        {
        // if(owner != rank){
          from_to_with_weight_t w;
          w.from = p;
          w.to = to;
          w.weight = weight;

          MPI_Isend(&w, 1, mpi_weight, owner, 0, MPI_COMM_WORLD, &requests[local_from]);
          send_to[owner]++;
        }
        // printf("END rank = %d owner = %d\n",rank, owner);
        // fflush(stdout);
        // MPI_Barrier(MPI_COMM_WORLD);

      }
      // MPI_Barrier(MPI_COMM_WORLD);
      // exit(1);
      // if(iter == 2) exit(1);


      MPI_Reduce_scatter(send_to, &receive_count, recounts, MPI_UINT32_T, MPI_SUM, MPI_COMM_WORLD);

      for(vertex_id_t i = 0; i < receive_count; i++) {
        from_to_with_weight_t w;
        MPI_Status status;
        MPI_Recv(&w,1, mpi_weight, MPI_ANY_SOURCE, 0,MPI_COMM_WORLD,&status);
        printf("!!!%d %u %u %lf\n", rank, w.from, w.to, w.weight);
        vertex_id_t local_form = w.from - first_vertex;
        if (full_edges[local_form].weight > w.weight){
          // printf("first_component_index = %u |  %u %u\n", VERTEX_TO_GLOBAL(i, G->n, G->nproc, G->rank), second_component_index, j);
          full_edges[local_form].to = w.to;
          full_edges[local_form].weight = w.weight;

          // printf("!!!first_component_index = %u |  %u %u\n", VERTEX_TO_GLOBAL(i, G->n, G->nproc, G->rank), second_component_index, j);
        }
      }
      // for(vertex_id_t i = 0; i < G->local_n; i++) {
      //
      //   vertex_id_t p = full_edges[find(full_edges,i,G)].parent;
      //   int owner = VERTEX_OWNER(p, G->n, G->nproc);
      //   send_to[owner]++;
      //   // printf("rank = %d f = %u |  %u\n", rank, 0, full_edges[i].edge);
      //
      //   // printf("rank = %d %u --- %u %u\n", rank, full_edges[find(full_edges,i,G)].parent,
      //   // find_in_edges(foreign_edges, full_edges[i].edge, G), full_edges[i].edge);
      //
      // }
      //



      // if(iter == 2) exit(1);
        for (vertex_id_t i=0; i < G->local_n; i++)
        {
            // printf("HERE 1 %d\n", omp_get_num_threads());
            // edge_id_t edge = full_edges[i].edge;
            vertex_id_t in_requests = 0, first_component_index, second_component_index;
            for(vertex_id_t i = 0; i < size; i++){
              all[i].from = UINT32_MAX;
              all[i].to = UINT32_MAX;
              recv_all[i].from = UINT32_MAX;
              recv_all[i].to = UINT32_MAX;
              all_result[i].from = UINT32_MAX;
              all_result[i].to = UINT32_MAX;
            }
            if (full_edges[i].weight != DBL_MAX && full_edges[i].parent == first_vertex + i)
            {
                printf("%d rank\n",rank );
                first_component_index = full_edges[find(full_edges, i,G)].parent;
                // vertex_id_t second_component_index = G->endV[edge]; //compute_second_component_index(full_edges, G, edge);

                second_component_index;
                vertex_id_t vertex =  full_edges[i].to; //G->endV[edge];
                if(vertex >= first_vertex && vertex < first_vertex + G->local_n )
                second_component_index =  full_edges[find(full_edges, G->endV[foreign_components.find(vertex)->second] - first_vertex, G)].parent;
                else
                second_component_index =   find_in_edges(foreign_edges, foreign_components.find(vertex)->second, G);
                if (first_component_index == second_component_index){
                  continue;
                }
                // changed = true;
                // for(vertex_id_t i = 0; i < size; i++){
                //   all[i].from = UINT32_MAX;
                //   all[i].to = UINT32_MAX;
                //   recv_all[i].from = UINT32_MAX;
                //   recv_all[i].to = UINT32_MAX;
                //   all_result[i].from = UINT32_MAX;
                //   all_result[i].to = UINT32_MAX;
                // }
                int owner = VERTEX_OWNER(second_component_index, G->n, G->nproc);
                all[owner].from = first_component_index;
                all[owner].to = second_component_index;
              }
                // memset(in_requests,0,G->local_n * sizeof(vertex_id_t));
                MPI_Alltoall(all, 1, mpi_edge, recv_all, 1, mpi_edge,MPI_COMM_WORLD);
              from_to_t f;
              f.from = f.to = UINT32_MAX;
              if (full_edges[i].weight != DBL_MAX && full_edges[i].parent == first_vertex + i)
              {
                for(vertex_id_t a = 0; a  < size;a++){
                  if(recv_all[a].to == first_component_index &&
                    !(recv_all[a].from == second_component_index && first_component_index < second_component_index)  ){
                    in_requests++;
                  }
                }


                if(in_requests == 0){
                  f.from = first_component_index;
                  f.to = second_component_index;
                  printf("%u request = %u  \n",  VERTEX_TO_GLOBAL(i, G->n, G->nproc, G->rank),
                  full_edges[i].to);
                }
              }
                MPI_Allgather(&f, 1, mpi_edge, all_result, 1, mpi_edge, MPI_COMM_WORLD);
              if (full_edges[i].weight != DBL_MAX && full_edges[i].parent == first_vertex + i)
              {

                for(vertex_id_t k=0;  k < size; k++){
                  vertex_id_t from = all_result[k].from;
                  vertex_id_t to = all_result[k].to;
                  // printf("PREFINISHED RANK = %d with k = %d %u -> %u \n", rank, k, from, to);

                  // if(from == UINT32_MAX) continue;
                  if(from != UINT32_MAX){

                    changed = true;
                  if( VERTEX_OWNER(from, G->n, G->nproc) == rank){
                      vertex_id_t local_from = from - first_vertex;
                      if( VERTEX_OWNER(to, G->n, G->nproc) == rank){
                        full_edges[local_from].parent = to;
                        full_edges[local_from].end = false;
                      }else{
                        auto to_pair = foreign_components.find(to);
                        edge_id_t edge = to_pair -> second;
                        full_edges[local_from].parent = to;
                        foreign_edges[edge].home_reference = from;
                      }
                      // printf("dawdwad %u\n", find(full_edges, from - first_vertex,G));
                  }else{
                    auto from_pair = foreign_components.find(from);
                    if(from_pair == foreign_components.end()) continue;
                    edge_id_t from_edge = from_pair -> second;

                    auto to_pair = foreign_components.find(to);
                    if(to_pair == foreign_components.end()){
                      foreign_components.insert(pair <vertex_id_t, edge_id_t>(to,from_edge));
                      G->endV[from_edge] = to;
                    }else{
                      foreign_edges[from_edge].parent = to_pair -> second;
                    }
                  }
                  }
                  // printf("FINISHED RANK = %d\n", rank);
                  // fflush(stdout);
                  // MPI_Barrier(MPI_COMM_WORLD);

                  // if(rank == 0){
                  //   printf("from = %u to = %u \n", all_result[k].from, all_result[k].to);
                  // }
                }


                // if(rank == 0)
                // printf("rank = %d ff= %u \n",rank, full_edges[0].parent);
                // fflush(stdout);
                // MPI_Barrier(MPI_COMM_WORLD);
                //
                // printf("rank = %d rr= %u \n",rank, find(full_edges, 0, G));
                // fflush(stdout);
                // MPI_Barrier(MPI_COMM_WORLD);
                // exit(0);

                // printf("!!!!! first_component_index = %u |  %u %u\n", VERTEX_TO_GLOBAL(i, G->n, G->nproc, G->rank), second_component_index, edge);

                // MPI_Send(&comparator, 1, mpi_edge, VERTEX_OWNER(second_component_index, G->n, G->nproc), MPI_UNITE_REQUEST, MPI_COMM_WORLD);
                // printf("SEND MAIN\n");
                // MPI_Recv(&receive_comparator, 1, mpi_edge, VERTEX_OWNER(second_component_index, G->n, G->nproc), MPI_UNITE_RESPONSE, MPI_COMM_WORLD, &status);
                // printf("RECV MAIN\n");

                // if()


                // if(unite_request(full_edges, first_component_index, second_component_index, G)){
                //   added_edges[edge] = true;
                // }
            }
        }
        // printf("HERWEW\n");
        // MPI_Barrier(MPI_COMM_WORLD);
        // printf("HERWEW2\n");


            // printf("HERE\n");

            // #pragma omp flush(continue_find)

        printf("FINAL %d\n", rank);

        // for(vertex_id_t i=0; i < G->local_n; i++){
        //  printf("%u request = %u edge = %u \n",  VERTEX_TO_GLOBAL(i, G->n, G->nproc, G->rank),   G->endV[full_edges[i].edge], full_edges[i].edge);
        // }

        // MPI_Barrier(MPI_COMM_WORLD);
        // exit(0);
        // for(vertex_id_t i=0; i < size; i++){
          // MPI_Alltoall(count_to,1, MPI_UINT32_T, receive_from,1, MPI_UINT32_T,MPI_COMM_WORLD);
          // vertex_id_t
          // for(vertex_id_t i=0; i < size; i++){
          //   printf("%d receive from %u = %u\n", G->rank,i, receive_from[i]);
          // }
          // /
          // / MPI_Isend(&count_to[i], 1, MPI_UINT32_T, i, 0, MPI_COMM_WORLD);
        // }

        // exit(0);
    //
    map <vertex_id_t,vertex_id_t> components_map;


    if(!changed){
      trees.clear();
      for(vertex_id_t i = 0; i < G->local_n;i++){
        vertex_id_t component = find(full_edges, i, G);
        if(i != component) continue;
        trees.push_back(vector<edge_id_t>());
        components_map.insert(pair <vertex_id_t, vertex_id_t>(component,trees.size() - 1));
      }
      for(edge_id_t i = 0; i < G->local_m; i++){
        edge_id_t edge = i;
        if(!added_edges[edge])
          continue;
        vertex_id_t component = find(full_edges, G->endV[edge], G);
        auto addr = components_map.find(component);
        trees[addr->second].push_back(edge_to_global(edge,G));
      }
      free(full_edges);
      free(added_edges);
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
