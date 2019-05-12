#include <vector>
#include <mpi.h>
#include <float.h>
#include <stdlib.h>
#include <assert.h>
#include "defs.h"
#define UINT32_MAX  (0xffffffff)
#define UINT64_MAX  (0xffffffffffffffff)

using namespace std;

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

vertex_id_t find(tree_vertex_t* tree_vertexes, vertex_id_t i)
{
  // if(VERTEX_OWNER(first_component_index, G->n, G->nproc)){
  //   i =
  // }
  vertex_id_t new_i;
  while(i != tree_vertexes[i].parent){
    new_i = tree_vertexes[i].parent;
    tree_vertexes[i].parent = tree_vertexes[new_i].parent;
    i = tree_vertexes[new_i].parent;
  }
  return i;
}

void assign_cheapest(full_edge_t* full_edges, graph_t*  G, vertex_id_t component_index,vertex_id_t second_component_index, edge_id_t j)
{
  if (full_edges[component_index].to == UINT32_MAX || G->weights[full_edges[component_index].edge] > G->weights[j]){
    full_edges[component_index].to = second_component_index;
    full_edges[component_index].edge = j;
    full_edges[component_index].proc = G->rank;
    full_edges[component_index].weight = G->weights[j];
  }
}

void unite(tree_vertex_t* tree_vertexes, vertex_id_t x, vertex_id_t y)
{
    // int xroot = find(subsets, x);
    // int yroot = find(subsets, y);
    // Attach smaller rank tree under root of high
    // rank tree (Union by Rank)
    if (tree_vertexes[x].rank < tree_vertexes[y].rank)
        tree_vertexes[x].parent = y;
    else if (tree_vertexes[x].rank > tree_vertexes[y].rank)
        tree_vertexes[y].parent = x;

    // If ranks are same, then make one as root and
    // increment its rank by one
    else
    {
        tree_vertexes[y].parent = x;
        tree_vertexes[x].rank++;
    }
}



extern "C" void* MST_boruvka(graph_t *G) {
    int rank = G->rank, size = G->nproc;
    vertex_id_t local_n = G->local_n, number_of_components = G->local_n, n = G->n;
    vertex_id_t* count_to = new vertex_id_t[size];
    // vertex_id_t  = new vertex_id_t[size];
    for(int i = 0; i < size; i++) {
      count_to[i] = 0;
    }
    for(vertex_id_t i = 0; i < G->local_n; i++) {
      for (edge_id_t j = G->rowsIndices[i]; j < G->rowsIndices[i+1]; j++) {
        count_to[VERTEX_OWNER(G->endV[j], G->n, G->nproc)]++;
      }
    }
    for(int i = 0; i < size; i++) {
      if(rank == 0)
        printf("count_to[%d]=%u\n", i, count_to[i]);
    }
    exit(0);

    // tree_vertex_t* tree_vertexes = (tree_vertex_t*)malloc(sizeof(tree_vertex_t) * G->n);
    // vector < edge_id_t > min_edges;
    // vector < vertex_id_t > min_edges_to;
    //
    // vertex_id_t* met_components = new vertex_id_t[n];
    // full_edge_t* full_edges = new full_edge_t[n];
    // full_edge_t* received_full_edges = new full_edge_t[n];

    // bool changed, changed_now;
    // for(vertex_id_t i=0; i < local_n;i++){
    //   // tree_vertexes[i].parent =  i;//VERTEX_TO_GLOBAL(i, G->n, G->nproc, G->rank);
    //   // tree_vertexes[i].rank = 0;
    //   // met_components[i] = UINT32_MAX;
    //
    // }


    // while(true){
    //   changed = false;
    //   for(vertex_id_t i=0; i < G->n;i++){
    //     full_edges[i].to = UINT32_MAX;
    //   }
    //
    //   for(vertex_id_t i = 0; i < G->local_n; i++) {
    //     for (edge_id_t j = G->rowsIndices[i]; j < G->rowsIndices[i+1]; j++) {
    //       // vertex_id_t first_component_index = find(tree_vertexes, VERTEX_TO_GLOBAL(i, G->n, G->nproc, G->rank));
    //       // vertex_id_t second_component_index = find(tree_vertexes, G->endV[j]);
    //       //
    //       // if (first_component_index == second_component_index){
    //       //   continue;
    //       // }
    //       //   assign_cheapest(full_edges, G, first_component_index, second_component_index, j);
    //       //   assign_cheapest(full_edges, G, second_component_index, first_component_index, j);
    //     }
    //   }
    //   int delta,from,to;
    //   MPI_Status status;
    //   for (vertex_id_t i=0; i < G->n; i++)
    //   {
    //       // Check if cheapest for current set exists
    //       if (full_edges[i].to != UINT32_MAX)
    //       {
    //           vertex_id_t first_component_index = find(tree_vertexes, i);
    //           vertex_id_t second_component_index = find(tree_vertexes, full_edges[i].to);
    //           if (first_component_index == second_component_index){
    //             continue;
    //           }
    //           changed = true;
    //
    //           unite(tree_vertexes, first_component_index, second_component_index);
    //           if(full_edges[i].proc == rank){
    //             min_edges.push_back(full_edges[i].edge);
    //             min_edges_to.push_back(G->endV[full_edges[i].edge]);
    //           }
    //
    //       }
    //   }
    //
    //   if(!changed){
    //     fflush(stdout);
    //     for(vertex_id_t i = 0; i < n; i++){
    //       vertex_id_t component = find(tree_vertexes, i);
    //       if(met_components[component] == UINT32_MAX){
    //         vertex_id_t size = trees.size();
    //         trees.push_back(vector<edge_id_t>());
    //         met_components[component] = size;
    //       }
    //     }
    //     for(edge_id_t i = 0; i < min_edges.size(); i++){
    //       vertex_id_t component = find(tree_vertexes, min_edges_to[i]);
    //       trees[met_components[component]].push_back( edge_to_global(min_edges[i],G));
    //
    //     }
    //     return &trees;
    //
    //   }
    // }



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
