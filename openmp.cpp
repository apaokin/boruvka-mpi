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

#include <map>

#define UINT32_MAX  (0xffffffff)
#define UINT64_MAX  (0xffffffffffffffff)
#define NUM_THREADS 4
#define NUM_THREADS_UNITE 4

using namespace std;
omp_nest_lock_t *vertex_locks;
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


int compare_and_swap(int* reg, int oldval, int newval)
{
  // ATOMIC();
  int old_reg_val = *reg;
  if (old_reg_val == oldval)
     *reg = newval;
  // END_ATOMIC();
  return old_reg_val;
}

// int update(vertex_id_t  )
//

vertex_id_t find(full_edge_t* full_edges, vertex_id_t i, graph_t* G)
{
  // if(VERTEX_OWNER(first_component_index, G->n, G->nproc)){
  //   i =
  // }
  // printf("rank = %d Want  find  i= %u\n", G->rank, i);
  // vertex_id_t new_i;
  // printf("find %u\n",i);
  // #pragma omp critical
  vertex_id_t new_i;
  {
  while(i != full_edges[i].parent){
          new_i = full_edges[i].parent;

          // #pragma omp critical
          // {
             __sync_val_compare_and_swap(&full_edges[i].parent, new_i,full_edges[new_i].parent);
            // if(full_edges[i].parent == new_i)
            //   full_edges[i].parent = full_edges[new_i].parent;
            // }
              i = full_edges[new_i].parent;
        // }
        // }
  }
  }
  // printf("found %u\n",i);

  return i;
}


void assign_cheapest(full_edge_t* full_edges, graph_t*  G, vertex_id_t i, edge_id_t j)
{
  // print_source_graph(G);
  // printf("\ni = %u, edge = %u, i2 = %u\n", i, j, G->endV[j]);

  // #pragma omp critical(cheap)
  {
    vertex_id_t first_component_index = find(full_edges, i, G);
    vertex_id_t second_component_index = find(full_edges, G->endV[j], G);
    //
    if (first_component_index != second_component_index){
      //   return;
      // }

      // if (find(full_edges, i, G) == find(full_edges, G->endV[j], G)){
      //   printf("WFADWAFWA\n");
      //   exit(1);;
      // }

      // printf("\nedge = %u, c1 = %u, c2 = %u\n", edge, first_component_index, second_component_index);
      // exit(1);
      // omp_set_nest_lock(&vertex_locks[first_component_index]);
      // __sync_synchronize(full_edges[first_component_index].edge);
      edge_id_t edge = full_edges[first_component_index].edge;

      // #pragma omp critical(cheap)
      // {
      //
      // }
      // printf("candidate edge %u for %u  with w = %lf; my e = %u my w = %lf\n", j, first_component_index, G->weights[j], edge,  G->weights[edge]);
      if (edge == UINT64_MAX || G->weights[edge] > G->weights[j]){
        // full_edges[first_component_index].vertex = second_vertex;
        full_edges[first_component_index].edge = j;
        // __sync_lock_test_and_set (&full_edges[first_component_index].edge, j);
        // printf("ASSIGNED edge %u for %u  with w = %lf\n", j, first_component_index, G->weights[j]);

        // full_edges[first_component_index].weight = G->weights[j];
        // // full_edges[second_vertex].vertex = first_vertex;
        // // full_edges[second_vertex].edge = j;
        // full_edges[second_vertex].weight = G->weights[j];

      }
      // omp_unset_nest_lock(&vertex_locks[first_component_index]);

    }
  }
}

inline bool unite_helper(full_edge_t* full_edges, vertex_id_t x, vertex_id_t y, omp_nest_lock_t* vertex_locks, vertex_id_t rank_x, vertex_id_t rank_y){
  // omp_set_nest_lock(&vertex_locks[x]);
  bool res, check;
  #pragma omp critical
  {
    // #pragma omp critical(unite)
    {
      // check = full_edges[x].parent != x || full_edges[y].parent != y;
      // full_edges[x].parent;
    }

    // printf("I WANT TO UNITE %u to %u from %d \n",x, y, omp_get_thread_num());
    // check = full_edges[x].parent != x || full_edges[x].rank != rank_x;
    vertex_id_t old_parent = full_edges[x].parent;
    vertex_id_t old_rank = full_edges[x].rank;
    // printf(old_rank)
    // if old.parent != v or old.rank != rank_v then return false
  // if(full_edges[x].parent != x || full_edges[y].parent != y){
  if(full_edges[x].parent != x || full_edges[x].rank != rank_x){
  // printf("1  %u-> %u \n",x, old_parent );
  // printf("2  %u-> %u \n",y, full_edges[y].parent );
  // #pragma omp critical
  // {
  //   check = (old_parent != x || old_rank != rank_x);
  // }
  // if(old_parent != x || old_rank != rank_x){

    // omp_unset_nest_lock(&vertex_locks[x]);
    // printf("FAILED TO UNITE %u to %u from %d \n",x, y, omp_get_thread_num());
    res = false;
  }else{

    // if (!__sync_bool_compare_and_swap(&full_edges[x].parent, old_parent,y)) {
    //   printf("EXIT\n");
    //   return false;
    // };
    full_edges[x].parent = y;
    // printf("%u\n", full_edges[x].parent);

    // full_edges[x].parent = y;
    if(full_edges[x].rank == full_edges[y].rank){
      // #pragma omp critical
      full_edges[y].rank++;
    }else{
      // #pragma omp critical
      // full_edges[x].rank = full_edges[y].rank;
    }
    // printf("SUCCEEDED TO UNITE %u to %u | %u %u rank_x=%u old=%u from %d \n",x, y, full_edges[x].parent, old_parent,  full_edges[x].rank, old_rank , omp_get_thread_num());
    // fflush(stdout);

    res  = true;
  }

  // fflush(stdout);
  }
  // omp_unset_nest_lock(&vertex_locks[x]);
  return res;
}


bool new_update(full_edge_t* full_edges, vertex_id_t v,vertex_id_t rank_v, vertex_id_t u, vertex_id_t rank_u){

  full_edge_t old = full_edges[v];
  if(old.parent != v || old.rank != rank_v)
   return false;

  full_edge_t new_struct;
  new_struct.parent = u;
  new_struct.rank = rank_u;// := { u, rank_u }

  return __sync_bool_compare_and_swap((long long*)&full_edges[v], *(long long*)&old, *(long long*)&new_struct);
}

bool unite(full_edge_t* full_edges, vertex_id_t x, vertex_id_t y, omp_nest_lock_t* vertex_locks)
{
    // printf("%u %u\n", x,y);
    // int xroot = find(subsets, x);
    // int yroot = find(subsets, y);
    // Attach smaller rank tree under root of high
    // rank tree (Union by Rank)
    // if(full_edges[x].parent != x || full_edges[y].parent != y)  return false;
    bool res;
    {

    if (full_edges[x].rank < full_edges[y].rank)
          res =  unite_helper(full_edges, x, y, vertex_locks, full_edges[x].rank, full_edges[y].rank);
          // res =  new_update(full_edges, x, full_edges[x].rank,  y, full_edges[x].rank);
    else if (full_edges[x].rank > full_edges[y].rank)
      // res =  new_update(full_edges, x, full_edges[x].rank,  y, full_edges[x].rank);

          res =  unite_helper(full_edges, y, x, vertex_locks, full_edges[x].rank, full_edges[y].rank);


    // If ranks are same, then make one as root and
    // increment its rank by one
    else{
    // if(res && full_edges[x].rank == full_edges[y].rank)
        // return new_update(full_edges, x, full_edges[x].rank,  x, full_edges[x].rank+1);

        res = unite_helper(full_edges, x, y, vertex_locks, full_edges[x].rank, full_edges[y].rank);

    }
    }
    return res;
    // return true;
}



void sort_vertexes(graph_t *G, vertex_id_t* foreign_vertexes, vertex_id_t &current_foreign_index ){
  vertex_id_t previous_vertex;
  vertex_id_t* buf_foreign_vertexes = new vertex_id_t[G->local_m];
  int rank = G->rank;

  memcpy(buf_foreign_vertexes, G->endV, sizeof(vertex_id_t)*G->local_m);
  // memcpy(buf_foreign_vertexes, G->endV, sizeof(vertex_id_t)*G->local_m);
  print_source_graph(G);
  sort(buf_foreign_vertexes, buf_foreign_vertexes + G->local_m);
  print_source_graph(G);

  previous_vertex = buf_foreign_vertexes[0];
  current_foreign_index = G->local_n;
  if(VERTEX_OWNER(previous_vertex, G->n, G->nproc) != rank){
    foreign_vertexes[current_foreign_index++] = buf_foreign_vertexes[0];
  }
  for(edge_id_t i=1; i < G->local_m; i++){
    if(previous_vertex != buf_foreign_vertexes[i]){
      if(VERTEX_OWNER(buf_foreign_vertexes[i], G->n, G->nproc) != rank){
        foreign_vertexes[current_foreign_index++] = buf_foreign_vertexes[i];
      }
      previous_vertex =  buf_foreign_vertexes[i];
    }
  }
  for(vertex_id_t i=0; i < G->local_n; i++){
    foreign_vertexes[i] = VERTEX_TO_GLOBAL(G->rowsIndices[i], G->n, G->nproc, G->rank);
  }

  vertex_id_t first_element = VERTEX_TO_GLOBAL(0, G->n, G->nproc, G->rank);
  for(vertex_id_t  i = 0; i < G->local_m;  i++){
    vertex_id_t  elem = G->endV[i];
    if(VERTEX_OWNER(elem, G->n, G->nproc) == rank){
      G->endV[i] = G->endV[i] - first_element;
    }
    else{
      vertex_id_t* index = lower_bound(&foreign_vertexes[G->local_n], &foreign_vertexes[current_foreign_index], elem);
      G->endV[i] = G->local_n + index - &foreign_vertexes[G->local_n];
      // if(rank == 3){
      //   printf("%u -- %u\n",elem,  G->endV[i]);
      // }
      // if(rank==0)
      //   printf("i=%u index=%u\n",i, index - foreign_vertexes);

    }
    // free(buf_foreign_vertexes);
  }

}

extern "C" void* MST_boruvka(graph_t *G) {
    int rank = G->rank, size = G->nproc;
    bool changed;
    vertex_locks = new omp_nest_lock_t[G->local_n];
    bool *added_edges = new bool[G->local_m];

    memset(added_edges,0,G->local_m * sizeof(bool));
    // for(vertex_id_t i = 0; i < G->local_n; i++) {
    //   for (edge_id_t j = G->rowsIndices[i]; j < G->rowsIndices[i+1]; j++) {
    //       vertex_id_t x = G->endV[j];
    //       if(i / 2 != x / 2)
    //         G->endV[j] = i;
    //   }
    // }
    // print_source_graph(G);

    // vertex_id_t local_n = G->local_n, number_of_components = G->local_n,
    // n = G->n, count_to = 0, foreign_vertexes_length = 0;
    // for(edge_id_t i=0; i < G->local_m; i++){

    full_edge_t *full_edges = new full_edge_t[G->local_n];
    // full_edge_t *full_edges = new full_edge_t[G->local_n];


    // int nitems = 4;
    // int blocklengths[nitems] = {1,1,1,1};
    // MPI_Datatype types[nitems] = {MPI_UINT32_T, MPI_UINT32_T, MPI_UINT64_T, MPI_DOUBLE};
    // MPI_Datatype mpi_edge;
    // MPI_Aint     offsets[nitems];
    //
    // offsets[0] = offsetof(full_edge_t, vertex);
    // offsets[1] = offsetof(full_edge_t, parent);
    // offsets[2] = offsetof(full_edge_t, edge);
    // offsets[3] = offsetof(full_edge_t, weight);
    //
    // MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_edge);
    // MPI_Type_commit(&mpi_edge);


    for(vertex_id_t i = 0; i < G->local_n; i++){
      full_edges[i].parent = i;
      full_edges[i].rank = 0;
    }

    for(vertex_id_t i = 0; i < G->local_n; i++){
      omp_init_nest_lock(&vertex_locks[i]);
    }
    // memcpy(buf_foreign_vertexes, G->endV, sizeof(vertex_id_t)*G->local_m);
    while(true){
        bool changed = false;
        for(vertex_id_t i = 0; i < G->local_n; i++){
          // full_edges[i].parent = i;
          full_edges[i].edge = UINT64_MAX;
          // full_edges[i].vertex = foreign_vertexes[i];
          // full_edges[i].rank = 0;

        }
        #pragma omp parallel for
        for(vertex_id_t i = 0; i < G->local_n; i++) {
          for (edge_id_t j = G->rowsIndices[i]; j < G->rowsIndices[i+1]; j++) {
              // buffer_to[begin_to[VERTEX_OWNER(G->endV[j], G->n, G->nproc)]] =
            // vertex_id_t first_component_index = find(full_edges,i, G);
            // vertex_id_t second_component_index = find(full_edges, G->endV[j], G);
            // // //
            // if (first_component_index == second_component_index){
            //   continue;
            // }
            {
              assign_cheapest(full_edges, G, i, j);
            }
          }
        }

      #pragma omp parallel for
      for (vertex_id_t i=0; i < G->local_n; i++)
      {
          // Check if cheapest for current set exists
          if (full_edges[i].edge != UINT64_MAX)
          {
              vertex_id_t first_component_index = find(full_edges, i,G);
              vertex_id_t second_component_index = find(full_edges, G->endV[full_edges[i].edge], G);
              if (first_component_index == second_component_index){
                continue;
              }
              changed = true;
              {
                // printf("I WANT TO UNITE %u and %u from %d \n",first_component_index, second_component_index, omp_get_thread_num());
                // fflush(stdout);
                if(unite(full_edges, first_component_index, second_component_index,vertex_locks)){
                  // printf("added edge %u\n",full_edges[i].edge);
                  // fflush(stdout);
                  // #pragma omp critical
                  added_edges[full_edges[i].edge] = true;
                  // {
                    // printf("THR %d\n", omp_get_thread_num());
                  // }
                  // printf("SUCCESS FINISHED %u and %u from %d \n",first_component_index, second_component_index, omp_get_thread_num());
                  // fflush(stdout);
                }
                else{
                  // printf("FAIL FINISHED %u and %u from %d \n",first_component_index, second_component_index,  omp_get_thread_num());
                  // fflush(stdout);

                }

              }
              // omp_unset_nest_lock(&vertex_locks[first_component_index]);
              // omp_unset_nest_lock(&vertex_locks[second_component_index]);


          }
      }
      // printf("ITER\n");
      // fflush(stdout);
    //
    map <vertex_id_t,vertex_id_t> components_map;


      if(!changed){
        // free(full_edges);
        // free(vertex_locks);

        // fflush(stdout);
        // for(vertex_id_t i = 0; i < G->local_n; i++){
        //   vertex_id_t component = find(tree_vertexes, i);
        //   if(met_components[component] == UINT32_MAX){
        //     vertex_id_t size = trees.size();
        //     trees.push_back(vector<edge_id_t>());
        //     met_components[component] = size;
        //   }
        // }
        trees.clear();
        // printf("trees.clear\n");
        for(vertex_id_t i = 0; i < G->local_n;i++){
          vertex_id_t component = find(full_edges, i, G);
          if(i != component) continue;
          trees.push_back(vector<edge_id_t>());
          components_map.insert(pair <vertex_id_t, vertex_id_t>(component,trees.size() - 1));
          // trees[trees.size() - 1].push_back(edge_to_global(edge,G));
        }
        // printf("COMP FINISHED\n");
        // fflush(stdout);

        // printf("INTERMED\n");
        // fflush(stdout);


        for(edge_id_t i = 0; i < G->local_m; i++){
          edge_id_t edge = i;
          if(!added_edges[edge])
            continue;
          // printf("BEFORE FIND %u max = %u \n", edge, UINT64_MAX);
          // fflush(stdout);
          // if(edge > 1000)
          // continue;
          vertex_id_t component = find(full_edges, G->endV[edge], G);
          // printf("edge %u component  %u\n", edge, component);
          auto addr = components_map.find(component);
          trees[addr->second].push_back(edge_to_global(edge,G));
          // for(edge_id_t i = 0; i < trees.size(); i++){
          //   printf("\ntree[%u]\n", i);
          //   for(edge_id_t j = 0; i < trees[i].size(); j++){
          //     printf("%u ", j);
          //   }
          // }

          // cout << trees << endl;
          //
          //

        }
        // printf("COMPLETED\n");
        // fflush(stdout);
        free(full_edges);
        free(vertex_locks);
        free(added_edges);
        // printf("FINISHED\n");
        // fflush(stdout);
        return &trees;
    //
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
