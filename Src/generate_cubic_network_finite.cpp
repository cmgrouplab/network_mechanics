//a random network with Possion distribution of nodes
//each node has 4 connection, to make sure mechanical stability


//author: yang.jiao.2@asu.edu
//modified: 07/02/2017
//this was originally wrote for Yongming to generate isostatic network, by perturbing
// a randomly diluted SC network, but the perturbation was not implemented correctly




using namespace std;


#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//configuration and basic properties


#define MAXL 10
#define num_node MAXL*MAXL*MAXL //number of node in the network
#define num_edge 16000 //number of edges in the network
#define max_neig 10 //the maximal number of neighbors

double coll_node[num_node][3]; //node configuration
int coll_edge[num_edge][2]; //edge configuration
int neig_ct[num_node]; //count the neighbor for each node for a neighbor list..
int neig_list[num_node][max_neig]; //list of neighbors, intialized when reading in the network

int connection_map[num_node][num_node]; //a connection map for easy handling of bonds
double dist_matrix[num_node][num_node]; //the dist map


double box_len = 1.0;
double cut_len = (box_len/(double)MAXL)*1.01;

ofstream fout;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/*
void init_data()
{
  coll
}
*/

/*
void generate_node()
{
  for(int i=0; i<num_node; i++)
    {
      for(int j=0; j<3; j++)
	coll_node[i][j] = box_len*(double)(rand()%100000)/(100000.0);
    }
}
*/


void generate_node_order()
{
  int temp_ct = 0;

  for(int i=0; i<MAXL; i++)
    for(int j=0; j<MAXL; j++)
      for(int k=0; k<MAXL; k++)
	{
	  coll_node[temp_ct][0] = (double)i*box_len/(double)MAXL;
	  coll_node[temp_ct][1] = (double)j*box_len/(double)MAXL;
	  coll_node[temp_ct][2] = (double)k*box_len/(double)MAXL;

	  temp_ct++;
        }
}

double dist(int index1, int index2)
{
  double temp_dist2 = 0.0;

  for(int i=0; i<3; i++)
    temp_dist2 += (coll_node[index1][i] - coll_node[index2][i])*(coll_node[index1][i] - coll_node[index2][i]);

  return sqrt(temp_dist2);
}

double dist_periodic(int index1, int index2)
{
  double temp_dist2 = 0.0;
  double temp_val = 0;

  for(int i=0; i<3; i++)
    {
       temp_val = fabs(coll_node[index1][i] - coll_node[index2][i]);   
       if(temp_val > 0.5*box_len) temp_val = box_len - temp_val;
       
       temp_dist2 += temp_val*temp_val;
    }
  return sqrt(temp_dist2);
}


void get_dist_matrix()
{
  for(int i=0; i<num_node; i++)
    for(int j=0; j<num_node; j++)
      dist_matrix[i][j] = dist(i, j);
}


void check_connection_order()
{
  get_dist_matrix();

  //init connection map
  for(int i=0; i<num_node; i++)
    for(int j=0; j<num_node; j++)
      connection_map[i][j] = 0;

  //reset of neighbor counters
  for(int i=0; i<num_node; i++)
    neig_ct[i] = 0;
  
  //check the dist matrix
  for(int i=0; i<num_node; i++)
    {
      cout<<"check node "<<i<<" for neighbors..."<<endl;

      for(int j=0; j<num_node; j++)
	{
	  if(dist_matrix[i][j]<cut_len && i!=j)
	    {
	      connection_map[i][j] = 1;

	      neig_ct[i]++;
	    }
	  
	}
    }

  
  for(int i=0; i<num_node; i++)
    {
      if(neig_ct[i]== 3)
	cout<<"neig "<<i<<" = "<<neig_ct[i]<<endl;
    }
  
}

//this was not implemneted correctly: the removed bonds are all along the z-direction
//the 3D connectivitiy was lost ..., but good as a testing case anyway...
void reduce_connection()
{
  for(int i=0; i<num_node; i++)
    {
      for(int j=0; j<num_node; j++)
	{
	  if(neig_ct[i]>4 && connection_map[i][j] == 1 && neig_ct[j]>4)
	    {
	      
	      if(connection_map[j][i] != 1)
		{
		  cout<<"connection is wrong"<<endl;
		  exit(1);
		}
	      
	      connection_map[i][j] = connection_map[j][i] = 0;
	      
	      neig_ct[i]--;
	      neig_ct[j]--;
	    }
	  
	}
    }

  //reset of neighbor counters
  for(int i=0; i<num_node; i++)
    neig_ct[i] = 0;
  
  //check the dist matrix
  for(int i=0; i<num_node; i++)
    {
      //cout<<"check node "<<i<<" for neighbors..."<<endl;

      for(int j=0; j<num_node; j++)
	{
	  if(connection_map[i][j] == 1)
	    neig_ct[i]++;
	}
    }

  for(int i=0; i<num_node; i++)
    {
      if(neig_ct[i]> 4 || neig_ct[i]==3)
	cout<<"neig "<<i<<" = "<<neig_ct[i]<<endl;
    }
}


void check_connection()
{
  get_dist_matrix();

  //init connection map
  for(int i=0; i<num_node; i++)
    for(int j=0; j<num_node; j++)
      connection_map[i][j] = 0;
  

  //reset of neighbor counters
  for(int i=0; i<num_node; i++)
    neig_ct[i] = 0;

  //check the dist matrix
  for(int i=0; i<num_node; i++)
    {
      cout<<"check node "<<i<<" for neighbors..."<<endl;

      for(int j=i+1; j<num_node; j++)
	{
	  if(neig_ct[i] == 4) break;

	  if(dist_matrix[i][j]<cut_len)
	    {
	      connection_map[i][j] = connection_map[j][i] = 1;

	      neig_ct[i]++;
	      neig_ct[j]++;
	    }
	  
	}

      if(neig_ct[i]<4)
	{
	  cout<<"cut_len too small! neig "<<i<<" = "<<neig_ct[i]<<endl;
	}
    }

  for(int i=0; i<num_node; i++)
    {
      if(neig_ct[i]!= 6)
	cout<<"neig "<<i<<" = "<<neig_ct[i]<<endl;
    }
    
}


void rand_move_node()
{
  double MOD = 1.5*cut_len;
  
  for(int i=0; i<num_node; i++)
    {
      for(int j=0; j<3; j++)
	{
	  coll_node[i][j] += MOD*((double)(rand()%100)/(100.0));
	  
	  /*
	  //make sure the periodic BC ....
	  if(coll_node[i][j]<0) coll_node[i][j] = coll_node[i][j] + box_len;
	  else if(coll_node[i][j]>box_len) coll_node[i][j] = coll_node[i][j] - box_len;
	  */
	}
    }
}

void print_node()
{
  cout<<"print node"<<endl;

  fout.open("polymer.nod");

  for(int i=0; i<num_node; i++)
    {
      fout<<(i+1)<<"\t"<<(coll_node[i][0])<<"\t"<<coll_node[i][1]<<"\t"<<coll_node[i][2]<<"\t"<<neig_ct[i]<<endl;
    }

  fout.close();
}

void print_edge()
{
  cout<<"print edge"<<endl;

  fout.open("polymer.edg");

  int temp_ct = 1;

  for(int i=0; i<num_node; i++)
    for(int j=(i+1); j<num_node; j++)
      {
	if(connection_map[i][j] == 1)
	  {
	    fout<<temp_ct<<"\t"<<(i+1)<<"\t"<<(j+1)<<endl;

	    temp_ct++;
	  }
      }

  fout.close();
}


void print_neigh()
{
     fout.open("polymer.lst");
     
     for(int i=0; i<num_node; i++)
     {
       fout<<(i+1)<<"\t\t";
         
         for(int j=0; j<num_node; j++)
            if(connection_map[i][j] == 1 && i!=j)
	           fout<<(j+1)<<"\t";
         
         fout<<endl;
     }
     
     fout.close();
}

main()
{
  srand(time(NULL));

  generate_node_order();

  check_connection_order();

  //reduce_connection();

  //rand_move_node();

  print_node();

  print_edge();
  
  print_neigh();
}
