//read in small periodic network configurations
//generate larger finite network configuration by replicating the smaller ones along different directions

//author: yang.jiao.2@asu.edu
//modified 02/24/16, two bugs found and corrected


//modified 06/28/2017
//mainly used to covert .edg file to .lst file
//also keeping the periodic BC

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

using namespace std;


#define N_node 5001
#define N_bond 8451
#define n_dim 3 //spatial dimension
#define n_max 30 //max number of neighbors
#define n_box 1 //number of box along each direction

double L = 1.0; //the periodicity of the length

#define N_node_F n_box*n_box*n_box*N_node
#define N_bond_F n_box*n_box*n_box*N_bond //since some bonds will be cut, the actual number will be smaller and a counter is required

double node_coords[N_node][n_dim];
int bond_list[N_bond][2]; //need to reduce the node index by 1 for connectivity check
int connect_map[N_node][N_node]; //this is the connectivity matrix, to map the nearest neighbor list
int neigh_num[N_node];
int neigh_list[N_node][n_max];


double new_node_coords[N_node_F][n_dim];
int new_bond_list[N_bond_F][2];
int new_neigh_num[N_node_F];
int new_neigh_list[N_node_F][n_max];

double dist0[N_node][n_max];

int bond_F_ct = 0; //the actual bond counter
int node_F_ct = 0;


ifstream fin;
ofstream fout;


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void read_config()
{
     //read in the nodes
     fin.open("coll.nod");
     double temp_x, temp_y, temp_z;
     int temp_num;
     
     for(int i=0; i<N_node; i++)
     {
        fin>>temp_num; //the index
        fin>>temp_x>>temp_y>>temp_z; //the coords
        fin>>temp_num; 
        
        node_coords[i][0] = temp_x;
        node_coords[i][1] = temp_y;
        node_coords[i][2] = temp_z;    
        
        neigh_num[i] = temp_num; //for the neighbor list 
     }
     
     fin.close();
     
     
     //now read in the connection, reduce the node index by 1, need to restore this when output the results
     //this is done by setting up a connection matrix
     for(int i=0; i<N_node; i++)
        for(int j=0; j<N_node; j++)
            connect_map[i][j] = 0;
            
     //now read in the connectivity 
     int index1, index2;
     
     fin.open("coll.edg");
     
     for(int i=0; i<N_bond; i++)
     {
         fin>>temp_num;
         fin>>index1>>index2;
         
         //reduce the index by one
         index1--;
         index2--;
         
         connect_map[index1][index2] = connect_map[index2][index1] = 1;
         
     }
     
     fin.close();
     
     //now get the neighbor list
     int temp_neigh_ct;
     
     for(int i=0; i<N_node; i++)
     {
         temp_neigh_ct = 0;
         
         for(int j=0; j<N_node; j++)
         {
             if(connect_map[i][j] == 1)
             {
                 
                 neigh_list[i][temp_neigh_ct] = j;
                 
                 temp_neigh_ct++;
             }
         }
         
         //another check
         if(temp_neigh_ct != neigh_num[i])
         {
             cout<<"neighbor number is wrong! node = "<<i<<endl;
             exit(1);                 
         }
     }
}

//do the translation to generate new coords and index
void generate_new_nodes()
{
     int temp_ct = 0;
     
     for(int i=0; i<n_box; i++)
       for(int j=0; j<n_box; j++)
         for(int k=0; k<n_box; k++)
            {
                 for(int m=0; m<N_node; m++)
                 {
                    new_node_coords[temp_ct][0] = node_coords[m][0] + k*L; 
                    new_node_coords[temp_ct][1] = node_coords[m][1] + j*L;
                    new_node_coords[temp_ct][2] = node_coords[m][2] + i*L;
                    
                    temp_ct ++;
                    
                    cout<<"generating node "<<temp_ct<<endl;
                 }
            }
            
     if(temp_ct != N_node_F)
     {
         cout<<"the total number of node is not correct!"<<endl;
         exit(1);
     }
     
}
 
 
double get_dist(int index1, int index2)
{
      double dx = new_node_coords[index1][0] -  new_node_coords[index2][0];  
      double dy = new_node_coords[index1][1] -  new_node_coords[index2][1]; 
      double dz = new_node_coords[index1][2] -  new_node_coords[index2][2];
      
      return sqrt(dx*dx + dy*dy + dz*dz);    
}
 
//check all images of neigh_index 
int check_neighbor(int node_index, int neigh_index)
{
    int min_neigh_index;
    double temp_min = 100000*L;
    double temp_dist;
    int temp_neigh_index;
    
    neigh_index = neigh_index%N_node; //tranlate this to the first box
    
    for(int i=0; i<(n_box*n_box*n_box); i++)
    {
            temp_neigh_index = neigh_index + i*N_node;
            temp_dist = get_dist(node_index, temp_neigh_index);
            
            if(temp_dist < temp_min)
            {
                temp_min = temp_dist;
                min_neigh_index = temp_neigh_index;
            }
    }
    
    if(temp_min < 0.25*L) return min_neigh_index;
    else return -1; 
} 
 
 
 
void generate_new_neigh_list()
{
     int temp_neigh_ct = 0;
     int temp_node_index;
     int temp_neigh_index;
     int new_neigh_index;
     
     int ave_neigh_num = 0;
     
     for(int i=0; i<N_node_F; i++)
     {
         temp_neigh_ct = 0;
         temp_node_index = i%N_node;
         
         for(int j=0; j<neigh_num[temp_node_index]; j++)  
         {
             temp_neigh_index = neigh_list[temp_node_index][j];
             
             new_neigh_index = check_neighbor(i, temp_neigh_index);
             
             if(new_neigh_index != -1)
             {
                 new_neigh_list[i][temp_neigh_ct] = new_neigh_index;
                 temp_neigh_ct ++;
             }
         }  
         
         new_neigh_num[i] = temp_neigh_ct;
         
         cout<<"neigh number of node "<<i<<" = "<<temp_neigh_ct<<endl;
         ave_neigh_num += temp_neigh_ct;
     }
     
     cout<<"ave_neigh_num = "<<(double)ave_neigh_num/((double)N_node_F)<<endl;

     //cout<<"ave_neigh_num ="<<ave_neigh_num<<endl;
} 
 

void print_node()
{
     fout.open("new_node.nod");
     
     for(int i=0; i<N_node_F; i++)
     {
       fout<<(i+1)<<"\t"<<new_node_coords[i][0]/((double)n_box*L)<<"\t"<<new_node_coords[i][1]/((double)n_box*L)<<"\t"<<new_node_coords[i][2]/((double)n_box*L)<<"\t"<<new_neigh_num[i]<<endl;
     }
     
     fout.close();
     
}

void print_neigh_list()
{
     fout.open("new_neigh.lst");
     
     for(int i=0; i<N_node_F; i++)
     {
       fout<<(i+1)<<"\t\t";
         
         for(int j=0; j<new_neigh_num[i]; j++)
	          fout<<(new_neigh_list[i][j]+1)<<"\t";
         
         fout<<endl;
     }
     
     fout.close();
}


void print_node_old()
{
     fout.open("new_node.nod");
     
     for(int i=0; i<N_node; i++)
     {
       fout<<(i+1)<<"\t"<<node_coords[i][0]<<"\t"<<node_coords[i][1]<<"\t"<<node_coords[i][2]<<"\t"<<neigh_num[i]<<endl;
     }
     
     fout.close();
     
}

void print_neigh_list_old()
{
     fout.open("new_neigh.lst");
     
     for(int i=0; i<N_node; i++)
     {
       fout<<(i+1)<<"\t\t";
         
         for(int j=0; j<neigh_num[i]; j++)
	         fout<<(neigh_list[i][j]+1)<<"\t";
         
         fout<<endl;
     }
     
     fout.close();
}


double get_dist_node(int index1, int index2)
{
     double temp_dist = 0.0;
     double temp_comp = 0.0;
     
     
     for(int i=0; i<n_dim; i++)
        {
           temp_comp = (node_coords[index1][i] - node_coords[index2][i]);
           if(fabs(temp_comp)> 0.5*L) temp_comp = L - fabs(temp_comp); 
           
           temp_dist += temp_comp*temp_comp;
        }
         
     return sqrt(temp_dist);
}

//get the iniital length of all fibers 
void get_dist0()
{
    int temp_index1, temp_index2;
    //this is done for each node, and its neighbors
    for(int i=0; i<N_node; i++)
    {
        temp_index1 = i;
        
        for(int j=0; j<neigh_num[temp_index1]; j++)
        {
                temp_index2 = neigh_list[temp_index1][j];
                
                dist0[i][j] = get_dist_node(temp_index1, temp_index2);
                
                if(dist0[i][j] > 0.5*L)
                  //cout<<"dist0 [ "<<i<<" "<<j<<endl;
                  cout<<"dist0 [ "<<i<<" "<<j<<" ] = "<<dist0[i][j]<<endl;
        }
    }
}


//do the index translation 
//void generate_new_bonds //when checking connections, need to check both nodes, to make sure 

//double check_connection(int index1, int index2)


//void get_connect_map()

//void print_large_network()

int main()
{

  read_config();

  //generate_new_nodes();

  //generate_new_neigh_list();

  print_node_old();

  print_neigh_list_old();

  get_dist0();
  //cout<<20000%N_bond<<endl;
    
  //int temp_num;
  //cin>>temp_num;

  return 0;
}

