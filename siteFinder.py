from ase import Atoms, Atom
from ase.io import read
from ase.geometry.analysis import Analysis
from ase.visualize import view
import numpy as np
import networkx as nx
import pandas as pd

def visualize_sites(sites, atoms):
    from ase import Atoms, Atom
    from ase.visualize import view
    for i in sites:
        atoms=atoms+Atoms('H',positions=[i])
    return view(atoms)
        
atoms=read('in.traj')

# Finding the exposed atoms

def top_layer_atoms(atoms):
    from ase.data import covalent_radii
    import numpy as np
    pos=atoms.positions
    slab_atom_types=atoms.get_atomic_numbers()
    top_layer=np.zeros(len(pos))==0
    for i in range(len(pos)):
        radii=covalent_radii[slab_atom_types[i]]
        for j in range(len(pos)):
            radij=covalent_radii[slab_atom_types[j]]
            if i!=j:
                D=pos[i]-pos[j]
                Dxy=(D[0]**2 + D[1]**2)**0.5
                if Dxy<radii+radij-0.4: #adjust the tolerance in xy plane
                    top_layer[i]=(top_layer[i]) and (pos[i][2]>=pos[j][2]-0.5) #adjust the tolerance in z direction
    return top_layer

def all_sites_possible(atoms):
    tops=atoms[top_layer_atoms(atoms)]
    pos=tops.get_positions()
    sites=[]
    #getting all possible adsorption sites: on top
    
    # On top
    
    for i in range(len(pos)):
        pos_copy=np.copy(pos[i])
        sites.append(pos_copy)
        sites[-1][2]=sites[-1][2]+2
    
    # Bridge
    
    ana=Analysis(tops)
    bonds=ana.all_bonds[0]
    for i in range(len(pos)):
        for j in bonds[i]:
            sites.append(((pos[i]+pos[j])/2))
            sites[-1][2]=sites[-1][2]+1.5  
            
    un=np.unique(np.array(sites),axis=0)
    return un

def spherical_neighbours(atoms,sites, neighbour_distance=None, n_neighbours=None):
    cell_param=np.array(np.sum(atoms.get_cell(),axis=0))
    cell_param[2]=0
    repited=tops.repeat((3,3,1))
    sites=sites+cell_param
    sites_relative_neighbours=[]
    if neighbour_distance is not None:
    # reletive spherical distance to the neighbors
        # chemical_symbs=np.unique(rep.get_chemical_symbols())
        for i in sites:
            for_this_site=[]
            for j in repited:
                chemsymb=j.symbol
                pos=j.position
                if np.linalg.norm(i-pos) <= neighbour_distance:
                    for_this_site.append([chemsymb, i-pos])
            sites_relative_neighbours.append(for_this_site)
        #now that we have the neighbours, we can convert them into spherical.
        for i in sites_relative_neighbours:
            for j in i:
                cart=j[1]
                r=np.linalg.norm(cart)
                theta=np.arctan(cart[1]/cart[0]) + (cart[0]<0)*np.pi
                if np.isnan(theta):
                    theta=np.arctan((cart[1]+0.0001)/cart[0]) + (cart[0]<0)*np.pi
                phi=np.arccos(cart[2]/r) + (cart[2]<0)*np.pi
                if np.isnan(phi):
                    phi=np.arccos((cart[2]+0.0001)/r) + (cart[2]<0)*np.pi
                j[1]=np.array([r,theta,phi])
        return sites_relative_neighbours
    # Specified number of neighbours
    elif n_neighbours is not None:
        pos=repited.positions
        chemsymb=np.array(repited.get_chemical_symbols())
        out=[]
        for i in sites:
            dist_from_this_site=np.linalg.norm(pos-i,axis=1)
            descending=np.argsort(dist_from_this_site)
            descending=descending[0:n_neighbours]
            neighbours=pos[descending]
            sites_relative_neighbours=i-neighbours
            
            symbs=chemsymb[descending]
            
            r=np.linalg.norm(sites_relative_neighbours-i,axis=1)
            theta=np.arctan((sites_relative_neighbours[:,1]+0.0001)/sites_relative_neighbours[:,0]) + (sites_relative_neighbours[:,0]<0)*np.pi
            phi=np.arccos((sites_relative_neighbours[:,2]+0.0001)/r) + (sites_relative_neighbours[:,2]<0)*np.pi
            coord=np.array([r,theta,phi]).T
            for o in range(len(symbs)):
                merged=[symbs[o],np.array(coord[o])]
            out.append(merged)

        return out

def sphere_to_cart(a):
    rho=a[:,0]
    sinphi=np.sin(a[:,2])
    cosphi=np.cos(a[:,2])
    sintheta=np.sin(a[:,1])
    costheta=np.cos(a[:,1])
    return np.array([rho*sinphi*costheta,rho*sintheta*sinphi,rho*cosphi]).T

# A function that we are going to use when checking the rotations.
# It checks if two sites are transformable to eachother with a rotation.
# == has not been used as a measure of equality of the coordinates. A tolerance is required
#  since a slight bit of difference (e.g. <0.05 Å or rad) is expected in most of the cases.

def compare_with_tolerance(check1, check2, T, tol=0.05,cart=True):
    if cart:
        ############take to cartesian###############
        #check1=sphere_to_cart(check1)
        check2=sphere_to_cart(check2)
        #T=sphere_to_cart(np.reshape(T,(1,3)))
        B_=check1+T
        B_=sphere_to_cart(B_)
        for i in B_:
            if not np.isin(1,np.sum(np.abs(check2-i)<tol,axis=1)//3):
                return False # means site1 and site2 are necessarly NOT equivalent
        return True # means site1 and site2 are transformable via T.
    else:
        # checks if all the rows in check1, can be transformed to sth in check2 by T.
        B_=check1+T
        for i in B_:
            if not np.isin(1,np.sum(np.abs(check2-i)<tol,axis=1)//3): #if one neighbour from site1 + T is not corresponding to one neighbour of site 2, then...
                return False # means site1 and site2 are necessarly NOT equivalent
        return True # means site1 and site2 are transformable via T.
    
# checking rotation and tossing away the duplicates
def detect_duplicates(neighbours):
    toss=[]
    for i in range(len(neighbors)):
        a=np.array(neighbors[i])
        types1=a.T[0]
        pos1=a.T[1]
        pos1=np.array([i for i in pos1])
        # find the neighbour with the lowest number of repetition to find the T
        min_occ=999999
        easiest_occ=''
        for occ in np.unique(types1):
            num_occ=np.count_nonzero(types1==occ)
            if num_occ<min_occ:
                min_occ=num_occ
                easiest_occ=occ
        # Now, the least frequent neighbour guides us toward finding the possible T.
        guide_from=pos1[types1==easiest_occ]
        for j in range(i+1,len(neighbors)):
            b=np.array(neighbors[j])
            types2=b.T[0]
            pos2=b.T[1]
            pos2=np.array([i for i in pos2])
            guide_to=pos2[types2==easiest_occ]
            for A in guide_from:
                for B in guide_to:
                    # A is a neighbour of the first site
                    # B is a neighbour of the second site
                    T=B-A
                    if abs(T[0])>0.05:
                        continue
                    similar=True
                    # if T[0]!=0:
                    #     continue
                    for occ in np.unique(types1):
                        check1=pos1[types1==occ]
                        check2=pos2[types2==occ]
                        # lets call a function that compares the two sites.
                        # check if A'+T is in {B'|B' is a neighbour of the second site.}
                        similar=similar and compare_with_tolerance(check1, check2, T, tol=0.7)
                        
                        if i==2 and j==27:# and occ=='O':
                            print('__________',i,j,'_______________')
                            print(occ)
                            print('T:',T)
                            print(compare_with_tolerance(check1, check2, T, tol=0.7))
                            print()
                        
                    if similar:
                        #print(i,'and',j,'are similar.\nThe T is:',T,'\n_____________________\n\n\n')
                        toss=np.append(toss,i)
                        break
                    #else:
                        #print(i,'and',j,'are NOT similar.','\n_____________________\n\n\n')
    return np.unique(toss)

                    
tops=atoms[top_layer_atoms(atoms)]
sites=all_sites_possible(tops)
# visualize_sites(sites,tops)

neighbors=spherical_neighbours(atoms,sites, neighbour_distance=7, n_neighbours=None)

tossables=detect_duplicates(neighbors)
tossables=[int(i) for i in tossables]
unique_sites=np.delete(sites,tossables,axis=0)
visualize_sites(unique_sites, atoms)
