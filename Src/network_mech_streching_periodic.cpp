//mechanics of collagen network
//using a force-based global relaxation, coupled with OpenMP 
//allow both linear and nonlinear collagen fiber model
//based on Long's code, rename the variables for better understanding


//started: 03/08/2016
//author: yang.jiao.2@asu.edu


//modified: 06/28/2017
//author: yang.jiao.2@asu.edu
//this original code is modified to model the stretching of semi-conducting fiber networks
//the cell is removed from the center, 
//uniaxial loading is applied, which is modeled by an affine transformation of the boundary 
//use periodic network configuration, more situable for this work




#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sstream>

using namespace std;

#define MAXX 50000 //make sure this number doesn't exceed the upper bound, set in rand()

#define n_dim 3 //spatial dimension
#define n_node 5001 //total number of nodes
#define n_max 25 //maxium number of neighbors

#define N_f 1000 //number of bins for force statsitics
#define N_r 30 //number of bins for decays of force/strain along radial direction

double Lx = 1.0; //linear size of the system along x direction
double Ly = 1.0; //linear size of the system along y direction
double Lz = 1.0; //linear size of the system along z direction
double L = (Lx + Ly + Lz)/3.0; //the characteristic size of the entire system

double bc_delta = 0.05*L; //the thickness of the boundary layer, in which the nodes are fixed, this is not useful for periodic BC

double coords[n_node][n_dim]; //coordinates of the nodes
int neigh_num[n_node]; //num of neighbors for each node
int neigh_list[n_node][n_max]; //the neighbor list

//int free_node_index[n_node]; //save the index of free neighbors
int free_node_flag[n_node]; //if free node = 1; otherwise = -1; this won't be used for periodic BC
//int n_free_node; //number of free nodes, will be saved in a separate array

double net_force[n_node][n_dim]; //the net force on each node
double force_residue[n_node]; //the magnitude of the net force
double ave_force_residue; //the average force residue, condition for converagence
double tran_mod = 0.001*L; //the translation mode, proportaional to the netforce on each node

double dist0[n_node][n_max]; //the original bond length associated with each node and its neighbors
double Ks = 1.0; //the springer constant
double gamma_buckle = 0.0; //the factor for buckling, the reduced modulus gamma_buckle*gamma
double lambda1 = 0.02; //the linear range of a fiber
double lambda2 = 0.03; //the exponential hardening length scale 

//these are left-over from the ECM-cell system
/*
double cell_rad = 0.2*L; //the radius of a spherical cell
//in generall, we should allow cell of any shape and morphology
double cell_contract_ratio = 0.5; //contraction ratio of a cell
double cell_center[n_dim]; //the center of the cell
*/

//parameters for relaxation
int N_loading = 1; //number for loading steps
double strain_xx = 0.0; // the normal strain
double strain_yy = 0.0;
double strain_zz = 0.0; 
//using different combination of the strain, to model different loading mode


int N_relax = 10; //number of relaxation after each loading

//statistics
double force_PDF[N_f]; //distribution of forces
double force_mod[n_node][n_max]; //magnitude of the forces along 
double force_max;
double force_min;

//for visiualizations 
#define N_vtk 150
int config[N_vtk][N_vtk][N_vtk]; //the configuration saveing the pixeled network
double force_threshold = 0.95; //to check agains f-f_min/f_max-f_min, to visualize the force chains


//for the force decay
double force_overall_decay[N_r]; //this is the averaged forced over all fibers
double force_chain_decay[N_r]; //this is the averaged over force-carrying fibers
int ct_overall[N_r];
int ct_chain[N_r];


//for the nonlinear strain
double lambda_decay[N_r]; //this is to characterize the nonlinear region
int ct_lambda[N_r];

//for input and output
ifstream fin;
ofstream fout;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//read in network and connections
void read_network()
{
     int temp_ind;
     
     fin.open("collagen.nod");
     for(int i=0; i<n_node; i++)
        {
             //read in the index, not very usefuall
             fin>>temp_ind;
             if(temp_ind != (i+1))
             {
                cout<<"index reading is not right for node coords!"<<endl;
                exit(1);
             }
             
             //now read in the coords
             for(int j=0; j<n_dim; j++)
                 fin>>coords[i][j];
                 
             //finally read in number of neighbors
             fin>>neigh_num[i];
        }     
     fin.close();
     
     //now read in the connections
     int self_neigh_flag = 0;
     int temp_ct = 0;
     fin.open("collagen.lst");
     for(int i=0; i<n_node; i++)
        {
            fin>>temp_ind;  
            if(temp_ind != (i+1))
             {
                cout<<"index reading is not right for neighbor list!"<<endl;
                exit(1);
             }

	    temp_ct = 0;           

            //now read in the neighbor index, make sure to substract 1
            for(int j=0; j<neigh_num[i]; j++)
               {
                    fin>>temp_ind;

                    if(temp_ind == (i+1))
		      {
                        cout<<"self-neighoring..."<<temp_ind<<endl;
                        self_neigh_flag = 1;
		      }
                    else
		      {
                        neigh_list[i][temp_ct] = temp_ind - 1;
                        temp_ct++;
		      }
               } 
             
            neigh_num[i] = temp_ct;
        }
     fin.close();
 }
 
void print_neigh(int index)
{
     cout<<"node "<<index<<" = ";
     for(int i=0; i<neigh_num[index]; i++)
         cout<<neigh_list[index][i]<<" ";
     cout<<endl;
        
 }
 

int check_boundary(int index)
{
    int temp_flag = 0;
    
    for(int i=0; i<n_dim; i++)
    {
       if(coords[index][i] < bc_delta || coords[index][i] > (L - bc_delta))
           temp_flag++;
    }
    
    if(temp_flag > 0) return -1;
    else return 1;
} 

 
//this is to be consistent with loading condition
//dim_index gives the direction of loading, and the boundaries to be fixed 
int check_boundary(int index, int dim_index) 
{
    int temp_flag = 0;
    double temp_L[n_dim];
    temp_L[0] = Lx;
    temp_L[1] = Ly;
    temp_L[2] = Lz;
    
    if(coords[index][dim_index] < bc_delta || coords[index][dim_index] > (temp_L[dim_index] - bc_delta))
           temp_flag++;
    
    if(temp_flag > 0) return -1;
    else return 1;
} 

 
 
//fix the boundary  
void get_boundary()
{
     //first, initialize the flag matrax 
     for(int i=0; i<n_node; i++)
        free_node_flag[i] = check_boundary(i, 0); //make sure the fixed boundary is consistent with the loading condition
}

 
//for the periodic boundary
void get_boundary_periodic()
{
     //in the case of periodic BC, all nodes are free nodes
     for(int i=0; i<n_node; i++)
        free_node_flag[i] = 1;
}
 
 
/*
double get_dist_cell(int index)
{
     double temp_dist = 0.0;
     
     for(int i=0; i<n_dim; i++)
         temp_dist += (coords[index][i] - cell_center[i])*(coords[index][i] - cell_center[i]);
         
     return sqrt(temp_dist);
} 
 */
 
/* 
//get the cell, can also read in complex shape and morphology 
//just set flag of noes within the cell to be -1 
void get_cell()
{
     //first, get the cell in the center
     for(int i=0; i<n_dim; i++)
         cell_center[i] = 0.5*L;
         
     //now check all nodes,
     int temp_ct = 0;
     for(int i=0; i<n_node; i++)
     {
         if(get_dist_cell(i) < cell_rad)
	   {
	     free_node_flag[i] = -1;
	     temp_ct++;
	   }
     }
     
     cout<<"cell_node = "<<temp_ct<<endl;
     
}
*/ 
 
/* 
//contract the node within the region specificed by radius, and also non-free nodes, i.e., with flag -1 
void get_contraction()
{
     double temp_dist;
     //now check all nodes, 
     for(int i=0; i<n_node; i++)
     {
         temp_dist = get_dist_cell(i);
             
         if(temp_dist < cell_rad && free_node_flag[i] == -1)
            {
                //now resacle the coords of the nodes
                for(int j=0; j<n_dim; j++)
                   coords[i][j] = cell_center[j] +  cell_contract_ratio*(coords[i][j] - cell_center[j]);
            }
     }
}
*/


//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//apply the deformation of the simulation box, specified by the strain components
void get_deformation()
{
     //frist update the linear edge lenght of the box
     Lx = Lx*(1 + strain_xx);
     Ly = Ly*(1 + strain_yy);
     Lz = Lz*(1 + strain_zz);
     
     //now update the coordinates of all the nodes
     for(int i=0; i<n_node; i++)
     {
         coords[i][0] = coords[i][0]*Lx;
         coords[i][1] = coords[i][1]*Ly;
         coords[i][2] = coords[i][2]*Lz;
     }

}

double f_abs(double val)
{
   if(val>=0) return val;
   else return -val;       
}

//&&&&&&&&&&&&&&&&&&&&&&&&&&
//compute the distance between two nodes
//periodic BC is used in this case, considering cuboid box with different edge length
double get_dist_node(int index1, int index2)
{
     double temp_dist = 0.0;
     double temp_comp = 0.0;
     double temp_L[n_dim];
     
     temp_L[0] = Lx;
     temp_L[1] = Ly;
     temp_L[2] = Lz;
     
     for(int i=0; i<n_dim; i++)
        {
           temp_comp = (coords[index1][i] - coords[index2][i]);
           if(f_abs(temp_comp)> 0.5*temp_L[i]) temp_comp = temp_L[i] - f_abs(temp_comp); 
           
           temp_dist += temp_comp*temp_comp;
        }
         
     return sqrt(temp_dist);
}

//get the iniital length of all fibers 
void get_dist0()
{
    int temp_index1, temp_index2;
    //this is done for each node, and its neighbors
    for(int i=0; i<n_node; i++)
    {
        temp_index1 = i;
        
        for(int j=0; j<neigh_num[temp_index1]; j++)
        {
                temp_index2 = neigh_list[temp_index1][j];
                
                dist0[i][j] = get_dist_node(temp_index1, temp_index2);
                
                /*
                if(dist0[i][j] > 0.05*L)
                  //cout<<"dist0 [ "<<i<<" "<<j<<endl;
                  cout<<"dist0 [ "<<temp_index1<<" "<<temp_index2<<" ] = "<<dist0[i][j]<<endl;
                  */
        }
    }
}


//compute the net force of each node
//compute the residue for each node
//compute the average residue for each node
//we can use OpenMP for this for very large networks 
void get_net_force_nonlinear() 
{
     double temp_force[n_dim]; //save the force due to one neighbor node
     double temp_force_mod; //the magntidue of the force between i and j
     double temp_net_force_mod; //the magnitude of the net force
     double temp_dist; //this is the distance between two nodes
     double temp_lambda; //this is the srain of the fiber
     int temp_index; //index of the neighbor nodes
     
     
     ave_force_residue = 0; //the total force residue
     
     //now loop over all nodes, and compute the forces between each neighbor pair
     for(int i=0; i<n_node; i++)
     {
         //make sure the old force is cleared
         for(int j=0; j<n_dim; j++)
            net_force[i][j] = 0;
         //for the magnitude
         temp_net_force_mod = 0;
            
         for(int j=0; j<neigh_num[i]; j++)
         {
             temp_index = neigh_list[i][j]; //get the neighbor index;
             
             temp_dist = get_dist_node(i, temp_index);
             temp_lambda = (temp_dist - dist0[i][j])/dist0[i][j];
                             
             //now compute the magnitude of the force
             if(temp_lambda <= 0) //this is compression buckliing
             {
                 temp_force_mod = 0.0;
             }
             else if(temp_lambda>0 && temp_lambda<lambda1) //this is linear elastic regime
             {
                  temp_force_mod = Ks*temp_lambda;
             }
             else if(temp_lambda>lambda1) //this is the strain hardening regime 
             {
                  temp_force_mod = Ks*lambda2*exp((temp_lambda-lambda1)/lambda2) + Ks*(lambda1-lambda2);
		  /*
		  if(temp_force_mod>1000000)
		    cout<<"temp_force_mod = "<<temp_force_mod<<" "<<i<<" "<<temp_index<<" temp_lambda = "<<temp_lambda<<"free flag "<<free_node_flag[i]<<" "<<free_node_flag[j]<<" dist ="<<temp_dist<<" dist0 = "<<dist0[i][j]<<" coords[index]"<<coords[temp_index][0]<<" "<<coords[temp_index][1]<<" "<<coords[temp_index][2]<<" coords[i]"<<coords[i][0]<<" "<<coords[i][1]<<" "<<coords[i][2]<<endl;
		  */
	     }
             
              //now we get the direction of the force, always pointing from i to neighbor, 
             for(int k=0; k<n_dim; k++)
                temp_force[k] = temp_force_mod*(coords[temp_index][k] - coords[i][k])/temp_dist;

	     //save the force mod
	     force_mod[i][j] = temp_force_mod;
                
             //assuming evertyhng is OK, add this force to the net force
             for(int k=0; k<n_dim; k++)
             {
                net_force[i][k] += temp_force[k];
             } 
         }
         
         for(int k=0; k<n_dim; k++)
            temp_net_force_mod += net_force[i][k]*net_force[i][k];
         
         force_residue[i] = sqrt(temp_net_force_mod);
         
         //now compute the force residue of the free nodes
         if(free_node_flag[i] != -1)
            ave_force_residue += force_residue[i];
         
     }
     
     ave_force_residue = ave_force_residue/(double)n_node;
     
}


//this is for the linear version
//void get_net_force_linear() 
//we can use OpenMP for this for very large network
void get_net_force_linear() 
{
     double temp_force[n_dim]; //save the force due to one neighbor node
     double temp_force_mod; //the magntidue of the force between i and j
     double temp_net_force_mod; //the magnitude of the net force
     double temp_dist; //this is the distance between two nodes
     double temp_lambda; //this is the srain of the fiber
     int temp_index; //index of the neighbor nodes
     
     
     ave_force_residue = 0; //the total force residue
     
     //now loop over all nodes, and compute the forces between each neighbor pair
     for(int i=0; i<n_node; i++)
     {
         //make sure the old force is cleared
         for(int j=0; j<n_dim; j++)
            net_force[i][j] = 0;
         //for the magnitude
         temp_net_force_mod = 0;
            
         for(int j=0; j<neigh_num[i]; j++)
         {
             temp_index = neigh_list[i][j]; //get the neighbor index;
             
             temp_dist = get_dist_node(i, temp_index);
             temp_lambda = (temp_dist - dist0[i][j])/dist0[i][j];
             //if(temp_lambda != 0)
	     //  cout<<"temp_lambda = "<<temp_lambda<<" dist0 = "<<dist0[i][j]<<" "<<i<<" "<<temp_index<<endl;
                             
             //now compute the magnitude of the force
             if(temp_lambda <= 0) //this is compression buckliing
             {
                 temp_force_mod = gamma_buckle*Ks*temp_lambda;
             }
             else if(temp_lambda>0) //this is linear elastic regime
             {
                  temp_force_mod = Ks*temp_lambda;
             }
           
              //now we get the direction of the force, always pointing from i to neighbor, 
             for(int k=0; k<n_dim; k++)
                temp_force[k] = temp_force_mod*(coords[temp_index][k] - coords[i][k])/temp_dist;

	     //save the temp_force_mod
	     force_mod[i][j] = temp_force_mod;
                
             //assuming evertyhng is OK, add this force to the net force
             for(int k=0; k<n_dim; k++)
             {
                net_force[i][k] += temp_force[k];
             } 
         }
         
         for(int k=0; k<n_dim; k++)
            temp_net_force_mod += net_force[i][k]*net_force[i][k];
         
         force_residue[i] = sqrt(temp_net_force_mod);
         //if(force_residue[i] != 0)
	 //   cout<<force_residue[i]<<endl;


         //now compute the force residue of the free nodes
         if(free_node_flag[i] != -1)
            ave_force_residue += force_residue[i];
         
     }
     
     ave_force_residue = ave_force_residue/(double)n_node;
     
}


//this is for the bending -- not sure how to implement yet ...
//void get_net_force_bending()

//apply the force-based displacement to relax network
//update each node
//we can use OpenMP for this 
void relax_network()
{
     //the displacement for each node
     double d_coords[n_dim];
     double d_mod = 0;
     
     //for periodic BC
     double temp_L[n_dim];
     temp_L[0] = Lx;
     temp_L[1] = Ly;
     temp_L[2] = Lz;
     
     //loop over all nodes, but only dispalce the free nodes
     for(int i=0; i<n_node; i++)
     {
         if(free_node_flag[i] != -1)
         {
	         d_mod = 0;
             for(int j=0; j<n_dim; j++)
	            {
		           d_coords[j] = tran_mod*net_force[i][j];

		           d_mod += d_coords[j]*d_coords[j];
		            //make sure everything is within the box...
		             /*
		              if(coords[i][j]>=L || coords[i][j]<0)
		              coords[i][j] = coords[i][j] - tran_mod*net_force[i][j];
		              */
	             }

             if(sqrt(d_mod)>2*tran_mod) //limits the maxium step size
	           {
		            for(int j=0; j<n_dim; j++)
                        coords[i][j] += 2*tran_mod*d_coords[j]/sqrt(d_mod);
	           }
	         else
	           {
		            for(int j=0; j<n_dim; j++)
		                coords[i][j] += d_coords[j];
	           }
	         
	         //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
             //now, need to make sure the displacement if consistent with periodic BC 
             for(int j=0; j<n_dim; j++)
                {
                     if(coords[i][j]<0) coords[i][j] = coords[i][j] + temp_L[j];
                     else if(coords[i][j] >= temp_L[j]) coords[i][j] = coords[i][j] - temp_L[j];
                 } 
	           
          }
     }
}


//get statistics of the forces
//the magnitude of forces along each fiber is in force_mod
void get_force_PDF()
{
  cout<<"computing force statistcs"<<endl;
  
  //loop over force_mod, to find out force_min and force_max
  force_min = 100000000.0;
  force_max = -100000000.0;

  for(int i=0; i<n_node; i++)
    {
      for(int j=0; j<neigh_num[i]; j++)
	{
	  if(free_node_flag[i] != -1)
	    {
	      if(force_mod[i][j]>force_max) force_max = force_mod[i][j];
	      if(force_mod[i][j]<force_min) force_min = force_mod[i][j];

	      //if(force_mod[i][j]>0) cout<<force_mod[i][j]<<endl;
	    }
	}
    }

  double force_bin = (force_max - force_min)/(double)N_f;

  //now get the statistics
  for(int i=0; i<N_f; i++)
    force_PDF[i] = 0.0;

  int temp_bin_num;
  int temp_ct = 0;

  for(int i=0; i<n_node; i++)
    {
      for(int j=0; j<neigh_num[i]; j++)
	{
	  if(free_node_flag[i] != -1)
	    {
	      temp_bin_num = (int)floor((double)(force_mod[i][j]-force_min)/(double)force_bin);
	      
	      if(temp_bin_num<N_f)
		{
		  //cout<<"force magnitude is not right, leading to wrong bin number"<<endl;
		  //exit(1);

		  force_PDF[temp_bin_num]++;
		  temp_ct++;
		}
	      
	 
	      
	    }
	}
    }
  
  //now print out the pdf
  fout.open("force_PDF.xls");
  for(int i=1; i<N_f; i++)
    {
      fout<<(force_min+force_bin*i)/force_max<<"\t"<<force_PDF[i]<<endl;
    }
  fout.close();

  fout.open("force_mod.xls");
  for(int i=0; i<n_node; i++)
    for(int j=0; j<neigh_num[i]; j++)
      if(free_node_flag[i] != -1)
	fout<<force_mod[i][j]<<endl;
  fout.close();

}


//print out the vtk for visualizing the force network
void print_force_network()
{
  int n_segment = 40;
  double l_pixel_x = (double)Lx/(double)N_vtk;
  double l_pixel_y = (double)Ly/(double)N_vtk;
  double l_pixel_z = (double)Lz/(double)N_vtk;
  
  double l_pixel = (double)L/(double)N_vtk;
  
  double temp_vect[n_dim];
  double temp_coords[n_dim];
  int index1, index2;
  
  double temp_L[n_dim];
  temp_L[0] = Lx;
  temp_L[1] = Ly;
  temp_L[2] = Lz;

  

  cout<<"initializing the vtk configuration, N_vtk = "<<N_vtk<<endl;
  //check each fiber, check each force, if greater than threshold, print out
  for(int i=0; i<N_vtk; i++)
    for(int j=0; j<N_vtk; j++)
      for(int k=0; k<N_vtk; k++)
	{
	  //put in the cell first
	  //temp_coords[0] = i*l_pixel;
	  //temp_coords[1] = j*l_pixel;
	  //temp_coords[2] = k*l_pixel;

      /*
	  double temp_dist = 0.0;
	  for(int m=0; m<n_dim; m++)
	    temp_dist += (temp_coords[m] - cell_center[m])*(temp_coords[m] - cell_center[m]);
	  temp_dist = sqrt(temp_dist);

	  if(temp_dist < cell_rad)
	    config[i][j][k] = 250;
	  else
	    config[i][j][k] = 0;
	    */
	    
	    config[i][j][k] = 0;
	}

  //loop over all nodes to check the bonds
  for(int i=0; i<n_node; i++)
    {
      //cout<<"generating coloar map, for node "<<i<<endl;

      index1 = i;

      for(int j=0; j<neigh_num[i]; j++)
	{
	  index2 = neigh_list[index1][j];

	  for(int k=0; k<n_dim; k++)
	  {
         /*     
         //need to correctly deal with peridoic bc here,
         double temp_vect_comp; 
         
         for(int a=-1; a<=1; a++)
            {
                 temp_vect_comp = coords[index2][k] - coords[index1][k] + a*temp_L[k];
                 
                 if(f_abs(temp_vect_comp)<0.5*temp_L[k])
                    {
                        temp_vect[k] = temp_vect_comp; 
                        break;
                    }
            }         
          */    
	    temp_vect[k] = coords[index2][k] - coords[index1][k];
      }
	  //now do the segement
	  //cout<<"generating the segment ..."<<endl;
	  
	  //this is for better and easy visualizaing
	  double temp_mod = 0; 
	      for(int k=0; k<n_dim; k++)
	         temp_mod += temp_vect[k]*temp_vect[k];        
	    
      //cout<<"dist "<<index1<<" "<<index2<<" = "<<sqrt(temp_mod)<<endl;
      
      
      //now, make sure peridoci connection is considered
      if(get_dist_node(index1, index2)>0.5*L) cout<<"dist "<<index1<<" "<<index2<<" = "<<sqrt(temp_mod)<<endl;
             
          if(sqrt(temp_mod)<0.25*L)
	  {
	  for(int m=0; m<n_segment; m++)
	    {
	      for(int k=0; k<n_dim; k++)
		     temp_coords[k] = coords[index1][k] + ((double)m/(double)n_segment)*temp_vect[k];

	      //digitized the fiber
	      int temp_x =(int)floor(temp_coords[0]/l_pixel_x);
	      int temp_y =(int)floor(temp_coords[1]/l_pixel_y);
	      int temp_z =(int)floor(temp_coords[2]/l_pixel_z);

	      double temp_force_ratio = (force_mod[i][j] - force_min)/(force_max - force_min);

           /*
	      //now check the magnitdue of the force
	      if(temp_force_ratio>force_threshold && free_node_flag[i] != -1)
            {
		       config[temp_x][temp_y][temp_z] = 255;
             }
	      else config[temp_x][temp_y][temp_z] = 100;
	      */
	      
	      //just visualizing the network, without worrying about the forces
	      //only print out network fibers 
          config[temp_x][temp_y][temp_z] = 100;
	      
	    }
     }
	    
	}
    }

  //now. print out the vtk
  cout<<"printing vtk file now ..."<<endl;
  fout.open("collagen_volume.vtk");

  fout<<"# vtk DataFile Version 2.0"<<endl;
  fout<<"Volume example"<<endl;
  fout<<"ASCII"<<endl;
  fout<<"DATASET STRUCTURED_POINTS"<<endl;
  fout<<"DIMENSIONS "<<N_vtk<<" "<<N_vtk<<" "<<N_vtk<<endl;
  fout<<"ASPECT_RATIO 1 1 1"<<endl; //<<Lx<<"\t"<<Ly<<"\t"<<Lz<<endl; //for small deformation, this should be fine
  fout<<"ORIGIN 0 0 0"<<endl;
  fout<<"POINT_DATA "<<(N_vtk)*(N_vtk)*(N_vtk)<<endl;
  fout<<"SCALARS volume_scalars UNSIGNED_CHAR 1"<<endl;
  fout<<"LOOKUP_TABLE default"<<endl;

  for(int i=0; i<N_vtk; i++)
    for(int j=0; j<N_vtk; j++)
      for(int k=0; k<N_vtk; k++)
	{
	  fout<<config[k][j][i]<<endl;
	}

  fout.close();

}

/*
void get_force_decay()
{
  double temp_dist;
  double dist_bin = (double)L/(2.0*N_r);
  int bin_num;
  int index1, index2;

  double temp_coords[n_dim];

  double force_threshold = 0.0001; //to check agains f-f_min/f_max-f_min
  
  cout<<"computing force decay ..."<<endl;

  //now clear the arrays
  for(int i=0; i<N_r; i++)
    {
      force_overall_decay[i] = 0;
      force_chain_decay[i] = 0;

      ct_overall[i] = 0;
      ct_chain[i] = 0;
    }

  //loop overall nodes and compute force statistics
  for(int i=0; i<n_node; i++)
    {
      index1 = i;

      //cout<<"node = "<<i<<endl;
      
      if(free_node_flag[i] != -1)
	{
	  for(int j=0; j<neigh_num[i]; j++)
	    {
	      index2 = neigh_list[i][j];
	      
	      temp_dist = 0.5*(get_dist_cell(index1) + get_dist_cell(index2));
	      bin_num = (int)floor(temp_dist/dist_bin);

	      double temp_force_ratio = (force_mod[i][j] - force_min)/(force_max - force_min);
	      
	      //add to the larger force
	      if(temp_force_ratio > force_threshold)
		{
		  if(bin_num < N_r)
		    {
		      force_chain_decay[bin_num] += force_mod[i][j];
		      ct_chain[bin_num]++;
		    }
		}

	      //add to overall forces
	      if(bin_num < N_r)
		{
		  force_overall_decay[bin_num] += force_mod[i][j];
		  ct_overall[bin_num]++;
		}
	      
	    }
	}
    }

  cout<<"print out the force decay profile"<<endl;

  double SA; 

  //now print out the statitics
  fout.open("force_overall_decay.txt");
  for(int i=1; i<N_r; i++)
    {
      SA = 4*3.1415926*(dist_bin*i)*(dist_bin*i)*dist_bin;
      
      if(ct_overall[i] == 0)
	fout<<dist_bin*i<<"\t"<<0<<endl;
      else
	fout<<dist_bin*i<<"\t"<<force_overall_decay[i]/((double)ct_overall[i]*SA)<<endl;
    }
  fout.close();

  fout.open("force_chain_decay.txt");
  for(int i=1; i<N_r; i++)
    {
      SA = 4*3.1415926*(dist_bin*i)*(dist_bin*i)*dist_bin;

      if(ct_chain[i] == 0)
	fout<<dist_bin*i<<"\t"<<0<<endl;
      else
	fout<<dist_bin*i<<"\t"<<force_chain_decay[i]/((double)ct_chain[i]*SA)<<endl;
    }
  fout.close();
  
}
*/


//get lambda decya

/*
void get_lambda_decay()
{
  double temp_dist;
  double dist_bin = (double)L/(2.0*N_r);
  int bin_num;
  int index1, index2;
  double temp_lambda;
  
  cout<<"computing lambda decay ..."<<endl;

  //now clear the arrays
  for(int i=0; i<N_r; i++)
    {
      lambda_decay[i] = 0;
      ct_lambda[i] = 0;
    }

  //loop overall nodes and compute force statistics
  for(int i=0; i<n_node; i++)
    {
      index1 = i;

      //cout<<"node = "<<i<<endl;
      
      if(free_node_flag[i] != -1)
	{
	  for(int j=0; j<neigh_num[i]; j++)
	    {
	      index2 = neigh_list[i][j];
	      
	      temp_dist = 0.5*(get_dist_cell(index1) + get_dist_cell(index2));
	      bin_num = (int)floor(temp_dist/dist_bin);

	      //now use the same temp_dist for displacement
	      temp_dist = get_dist_node(i, index2);
	      temp_lambda = (temp_dist - dist0[i][j])/dist0[i][j];
                             
	      
	      //add to the larger force
	      if(temp_lambda > lambda1)
		{
		  if(bin_num < N_r)
		    {
		      lambda_decay[bin_num] += temp_lambda;
		      ct_lambda[bin_num]++;
		    }
		}
	      
	    }
	}
    }

  cout<<"print out the lambda decay profile"<<endl;

  double SA; 

  //now print out the statitics
  fout.open("lambda_decay.txt");
  for(int i=1; i<N_r; i++)
    {
      SA = 4*3.1415926*(dist_bin*i)*(dist_bin*i)*dist_bin;
      
      if(ct_lambda[i] == 0)
	fout<<dist_bin*i<<"\t"<<0<<endl;
      else
	fout<<dist_bin*i<<"\t"<<lambda_decay[i]/((double)ct_lambda[i]*SA)<<endl;
    }
  fout.close();
  
}
*/

//also need to print the flag?
//separate each deformation iteration is necessary
void print_network(int iteration)
{
     char fn [100];
     snprintf (fn, sizeof fn, "new_node_%02d.nod", (iteration+1));
     fout.open(fn);
     
     //print out the coords
     //fout.open("new_node.nod");
     for(int i=0; i<n_node; i++)
     {
       fout<<(i+1)<<"\t"<<coords[i][0]<<"\t"<<coords[i][1]<<"\t"<<coords[i][2]<<"\t"<<neigh_num[i]<<endl;
     }
     fout.close();
     
     
     //~~~~~~~~~~~~~~~~~~~~~~~~~~
     snprintf (fn, sizeof fn, "new_neigh_%02d.lst", (iteration+1));
     fout.open(fn);
     
     //print out the neighors
     //fout.open("new_neigh.lst");
     int temp_ind;
     for(int i=0; i<n_node; i++)
     {
       fout<<(i+1)<<"\t\t";   
       for(int j=0; j<neigh_num[i]; j++)
       { 
         //&&&&&&&&&&&&&&& a bug was found and corrected....
         temp_ind = neigh_list[i][j];
	     fout<<(temp_ind+1)<<"\t";
       }
       fout<<endl;
     }
     fout.close();
}




int main()
{
    //initialization
    read_network(); //read in the nodes and connections
    
    //print_neigh(0);
    
    
    get_dist0(); // get the original length of the fibers
    
    //get_boundary_periodic(); //fix the box boundary
    get_boundary();
    
    
    //print_neigh(0);
    
   // get_cell(); //put in the cell
    
    //now start iteration 
    for(int i=0; i<N_loading; i++)
    {
        get_deformation(); //contract the cell by the given ratio
        //print_neigh(0);
     
        
        //get_net_force_nonlinear(); //compute the net force
	    get_net_force_linear(); //compute the net force
        //print_neigh(0);
        
        for(int j=0; j<N_relax; j++)
        {
            relax_network();
            //print_neigh(0);
            
            //get_net_force_nonlinear(); //compute the net force
	        get_net_force_linear(); //compute the net force
            //print_neigh(0);
            
            
            cout<<"Deformation "<<i<<" Relax "<<j<<" ave_force_residue = "<<ave_force_residue<<endl;
        }
        
       print_network(i);
    }
    
    
    
    //print the force network, to check
    cout<<"printing vtk for force network..."<<endl;
    print_force_network();
    
    //now print out the staitics
    get_force_PDF();
   

    //get_force_decay();

    //get_lambda_decay();
    
    return 0;
}
