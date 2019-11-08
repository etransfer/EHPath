print(' ')
print('*****************************************')
print('Welcome! EHPath is able to analyze and rank electron/hole hopping pathways between an electron/hole donor and an electron/hole acceptor according to the overall mean residence time.')
print('*****************************************')
print(' ')

# get input file names
donor_name=input('Please enter the name of donor file (eg: donor.csv): ')
bridge_name=input('Please enter the name of bridge file, (eg: bridge.csv): ')
acceptor_name=input('Please enter the name of acceptor file, (eg: acceptor.csv): ')

donor_num=int(input('Please enter residue number of donor: '))
cutoff_num=int(input('Please enter cutoff_num (cutoff_num + 1 = maximum number of nodes in the pathway, Warning: A larger cutoff_num will cost more memory and leads to longer computation time): '))
total_paths=int(input('Please specify the number of pathways to be printed: '))
type_da=input('Please enter the type of donor/acceptor (electron or hole): ')
alpha_reorg=int(input('Please specify the value for α: '))
print(' ')
print('Start calculation........')

import pandas as pd
import networkx as nx

acceptors = pd.read_csv(acceptor_name,names=['atom','atomnumber','atomtype','residuetype','chain','residue or node number','x','y','z','a','b','c'])
donor = pd.read_csv(donor_name,names=['atom','atomnumber','atomtype','residuetype','chain','residue or node number','x','y','z','a','b','c'])
bridges = pd.read_csv(bridge_name,names=['atom','atomnumber','atomtype','residuetype','chain','residue or node number','x','y','z','a','b','c'])

##############Section for truncating redox groups##############
###Charge Acceptors###
#Truncated Tyr, Met, and Trp are represented as phenol, thioether, and indole respectively.
acceptors = acceptors[~(acceptors['atomtype'].isin(['N','CA','C','O','CB']) & (acceptors['residuetype'].isin(['TYR','MET','TRP'])) )]
#Truncated Cys is represented as methanethiol.
acceptors = acceptors[~(acceptors['atomtype'].isin(['N','CA','C','O']) & (acceptors['residuetype'].isin(['CYS'])) )]
#Truncated guanine and adenine, where only the nucleobase without the DNA backbone is considered. 
acceptors = acceptors[~(acceptors['atomtype'].isin(['P','OP1','OP2','O5\'','C5\'','C4\'','O4\'','C3\'','O3\'','C2\'','O2\'','C1\'']) & (acceptors['residuetype'].isin(['G','DG','DA','A'])) )]
acceptors = acceptors[~(acceptors['atomtype'].isin(['PG','O1G','O2G','O3G','O3B','PB','O1B','O2B','O3A','PA','O1A','O2A','O5\'','C5\'','C4\'','O4\'','C3\'','O3\'','C2\'','O2\'','C1\'']) & (acceptors['residuetype'].isin(['GTP'])) )]
#Heme
acceptors = acceptors[~(acceptors['atomtype'].isin(['CBB','CAB','CMB','CMA','CAA','CBA','CGA','O1A','O2A','O2D',
                                       'O1D','CGD','CBD','CAD','CMD','CAC','CBC','CMC']) & (acceptors['residuetype'].isin(['HEM'])) )]
###Charge Donors###
#Heme hole donor (cytochrome c peroxidase and cytochrome p450).
donor = donor[~(donor['atomtype'].isin(['CBB','CAB','CMB','CMA','CAA','CBA','CGA','O1A','O2A','O2D',
                                       'O1D','CGD','CBD','CAD','CMD','CAC','CBC','CMC']) & (donor['residuetype'].isin(['HEM'])) )]
#Gly hole donor (for BSS), where only the CA atom is considered.
donor = donor[~(donor['atomtype'].isin(['N','C','O']) & (donor['residuetype'].isin(['GLY'])) )]
#Truncated Tyr, Met, Trp, Cys, guanine, and adenine.
donor = donor[~(donor['atomtype'].isin(['N','CA','C','O','CB']) & (donor['residuetype'].isin(['TYR','MET','TRP'])) )]
donor = donor[~(donor['atomtype'].isin(['N','CA','C','O']) & (donor['residuetype'].isin(['CYS'])) )]
donor = donor[~(donor['atomtype'].isin(['P','OP1','OP2','O5\'','C5\'','C4\'','O4\'','C3\'','O3\'','C2\'','O2\'','C1\'']) & (donor['residuetype'].isin(['G','DG','DA','A'])) )]
donor = donor[~(donor['atomtype'].isin(['PG','O1G','O2G','O3G','O3B','PB','O1B','O2B','O3A','PA','O1A','O2A','O5\'','C5\'','C4\'','O4\'','C3\'','O3\'','C2\'','O2\'','C1\'']) & (donor['residuetype'].isin(['GTP'])) )]
###Bridging Sites###
#Truncated Tyr, Met, Trp
bridges = bridges[~(bridges['atomtype'].isin(['N','CA','C','O','CB']) & (bridges['residuetype'].isin(['TYR','MET','TRP'])) )]
#Truncated Cys
bridges = bridges[~(bridges['atomtype'].isin(['N','CA','C','O']) & (bridges['residuetype'].isin(['CYS'])) )]
#Heme (can be applied to studying multi-heme proteins, for example, MtrF)
bridges = bridges[~(bridges['atomtype'].isin(['CBB','CAB','CMB','CMA','CAA','CBA','CGA','O1A','O2A','O2D',
                                       'O1D','CGD','CBD','CAD','CMD','CAC','CBC','CMC']) & (bridges['residuetype'].isin(['HEM'])) )]

##############Section for defining redox groups and their associated charge transfer parameters##############
#self.Lamb is equal to one-half of the self-exchange reorganization energies, self.E is the oxidation potential, while self.r is the effective radius.
class node(object):
    
    def __init__(self, nodeNum, Coords,Residuetype,AtomQuant):
        self.nodeNum=nodeNum
        self.Coords=Coords
        self.Residuetype=Residuetype
        if self.Residuetype in ['GTP','G','DG']:
            self.Lamb=0.373
            self.E=1.29
            self.r=2.10 #radius
        if self.Residuetype in ['DA','A']:
            self.Lamb=0.2115
            self.E=1.42
            self.r=2.22 
        if self.Residuetype=='CYS':
            self.Lamb=1.27
            self.E=0.92
            self.r=1 #fictitious radius
        if self.Residuetype=='MET':
            self.Lamb=1.08
            self.E=1.66
            self.r=1
        if self.Residuetype=='TYR':
            self.Lamb=1.02
            self.E=0.93
            self.r=1
        if self.Residuetype=='TRP':
            self.Lamb=0.95
            self.E=1.02
            self.r=1
        if self.Residuetype=='SF4':
            self.Lamb=0.75
            self.E=0.08
            self.r=1
        if self.Residuetype=='HEM':
            self.Lamb=0.85
            self.E = 1.00
            self.r = 1
        if self.Residuetype=='GLY':
            self.Lamb=1.00
            self.E = 1.9
            self.r = 1
        self.AtomQuant=AtomQuant

##############Section for constructing the directed graph and the associated nodes##############
G=nx.DiGraph()

donor_list = list(donor['residue or node number'].unique())
acceptor_list=list(acceptors['residue or node number'].unique())
bridge_list=list(bridges['residue or node number'].unique())
nodes=[]

for e in donor_list:
    nodes.append(node(e,
                      donor[donor['residue or node number']==e][['x','y','z']],
                      max(donor[donor['residue or node number']==e]['residuetype']),
                      len(donor[donor['residue or node number']==e]),                      
                           ))

for e in bridge_list:
    nodes.append(node(e,
                      bridges[bridges['residue or node number']==e][['x','y','z']],
                      max(bridges[bridges['residue or node number']==e]['residuetype']),
                      len(bridges[bridges['residue or node number']==e]),                      
                           ))

for e in acceptor_list:
    nodes.append(node(e,
                      acceptors[acceptors['residue or node number']==e][['x','y','z']],
                      max(acceptors[acceptors['residue or node number']==e]['residuetype']),
                      len(acceptors[acceptors['residue or node number']==e]),                      
                           ))
for e in nodes:
    G.add_node(e.nodeNum)

##############Calculation of geometric center-to-center distance and shortest distance between two redox groups##############
def calculate_distance(a,b):
    return ((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2)**0.5
    
def get_distances(coordinate_a,coordinate_b):
    ax=list(coordinate_a['x'])
    ay=list(coordinate_a['y'])
    az=list(coordinate_a['z'])
    acoord=[(ax[i],ay[i],az[i]) for i in range(len(ax))]
    
    bx=list(coordinate_b['x'])
    by=list(coordinate_b['y'])
    bz=list(coordinate_b['z'])
    bcoord=[(bx[i],by[i],bz[i]) for i in range(len(bx))]
    
    mean_distance=calculate_distance((np.mean(ax),np.mean(ay),np.mean(az)),(np.mean(bx),np.mean(by),np.mean(bz)))
    
    all_combine=[[x,y] for x in acoord for y in bcoord]
    all_distance=[]
    for e in all_combine:
        all_distance.append(calculate_distance(e[0],e[1]))
    shortest_distance=min(all_distance)
    return mean_distance,shortest_distance

import numpy as np
def smallK(source,target): 
    source_node=[e for e in nodes if e.nodeNum==source][0]
    target_node=[e for e in nodes if e.nodeNum==target][0]
    mean_distance,shortest_distance=get_distances(source_node.Coords,target_node.Coords)

##############Section for defining interactions (i.e. Marcus parameters and rates) between two nodes (i.e. redox groups) in any hopping pathway##############        
##############Donor-Bridge##############
    
    #Compute Marcus parameters between the hole donor node (heme) and a bridge node (amino acid) for the forward HT step.  
    if type_da=='hole' and (source_node.nodeNum in donor_list and source_node.Residuetype in ['HEM']) and (target_node.nodeNum in bridge_list and target_node.Residuetype in ['CYS','MET','TRP','TYR']):
        Lambda=(source_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-17.45)/(mean_distance*17.45)/1.889726)+(target_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-12.8)/(mean_distance*12.8)/1.889726)
        Lambda=0.8+alpha_reorg*(Lambda*27.2114-0.8)
        deltaG=target_node.E-source_node.E
        VBB=2.7/(target_node.AtomQuant*target_node.AtomQuant)**0.5 * np.exp(-0.72*shortest_distance)
        VDD=np.mean(np.array([8.0])*np.exp(-1*np.array([2.6])*((shortest_distance + 2*3.63)-16.5)/2))*0.00012398
        VIF=(VDD*VBB)**0.5

    #Compute Marcus parameters between the hole donor node (heme) and a bridge node (amino acid) for the backward HT step.      
    if type_da=='hole' and (source_node.nodeNum in bridge_list and source_node.Residuetype in ['CYS','MET','TRP','TYR']) and (target_node.nodeNum in donor_list and target_node.Residuetype in ['HEM']):
        Lambda=(target_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-17.45)/(mean_distance*17.45)/1.889726)+(source_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-12.8)/(mean_distance*12.8)/1.889726)
        Lambda=0.8+alpha_reorg*(Lambda*27.2114-0.8)
        deltaG=target_node.E-source_node.E
        VBB=2.7/(source_node.AtomQuant*source_node.AtomQuant)**0.5 * np.exp(-0.72*shortest_distance)
        VDD=np.mean(np.array([8.0])*np.exp(-1*np.array([2.6])*((shortest_distance + 2*3.63)-16.5)/2))*0.00012398
        VIF=(VDD*VBB)**0.5             

    #Compute Marcus parameters between the hole donor node (heme) and a bridge node (heme) for the forward and backward HT step.  
    if type_da=='hole' and ((source_node.nodeNum in donor_list and source_node.Residuetype in ['HEM']) and (target_node.nodeNum in bridge_list and target_node.Residuetype in ['HEM'])) or ((source_node.nodeNum in bridge_list and source_node.Residuetype in ['HEM']) and (target_node.nodeNum in donor_list and target_node.Residuetype in ['HEM'])):
        Lambda=(source_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-17.45)/(mean_distance*17.45)/1.889726)+(target_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-17.45)/(mean_distance*17.45)/1.889726)
        Lambda=0.8+alpha_reorg*(Lambda*27.2114-0.8)
        deltaG=target_node.E-source_node.E
        VBB=np.mean(np.array([8.0])*np.exp(-1*np.array([2.6])*((shortest_distance + 2*3.63)-16.5)/2))*0.00012398
        VDD=np.mean(np.array([8.0])*np.exp(-1*np.array([2.6])*((shortest_distance + 2*3.63)-16.5)/2))*0.00012398
        VIF=(VDD*VBB)**0.5

    #For BSS: Compute Marcus parameters between the hole donor node (Gly) and a bridge node (amino acid) for the forward and backward HT step.  
    if type_da=='hole' and ((source_node.nodeNum in donor_list and source_node.Residuetype in ['GLY']) and (target_node.nodeNum in bridge_list and target_node.Residuetype in ['CYS','MET','TRP','TYR'])) or ((source_node.nodeNum in bridge_list and source_node.Residuetype in ['CYS','MET','TRP','TYR']) and (target_node.nodeNum in donor_list and target_node.Residuetype in ['GLY'])):
        Lambda=(source_node.Lamb+target_node.Lamb)/27.2114+(1/2.2-1/4)*(mean_distance-12.8)/(mean_distance*12.8)/1.889726
        Lambda=0.8+alpha_reorg*(Lambda*27.2114-0.8)
        deltaG=target_node.E-source_node.E
        VIF=2.7/(source_node.AtomQuant*target_node.AtomQuant)**0.5 * np.exp(-0.72*shortest_distance)
 
    #For primase-DNA: Compute Marcus parameters between the *electron* donor node ([4Fe4S]) and a bridge node (amino acid) for the forward *ET* step.     
    if type_da=='electron' and (source_node.nodeNum in donor_list and source_node.Residuetype in ['SF4']) and (target_node.nodeNum in bridge_list and target_node.Residuetype in ['CYS','MET','TRP','TYR']):
        Lambda=(source_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-14.44)/(mean_distance*14.44)/1.889726)+(target_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-12.8)/(mean_distance*12.8)/1.889726)
        Lambda=0.8+alpha_reorg*(Lambda*27.2114-0.8)
        deltaG=source_node.E-target_node.E
        VBB=2.7/(target_node.AtomQuant*target_node.AtomQuant)**0.5 * np.exp(-0.72*shortest_distance)
        VDD=10**(1.73-0.42*(shortest_distance+2*1.87))
        VIF=(VDD*VBB)**0.5
        
    #For primase-DNA: Compute Marcus parameters between the *electron* donor node ([4Fe4S]) and a bridge node (amino acid) for the backward *ET* step.     
    if type_da=='electron' and (source_node.nodeNum in bridge_list and source_node.Residuetype in ['CYS','MET','TRP','TYR']) and (target_node.nodeNum in donor_list and target_node.Residuetype in ['SF4']):
        Lambda=(target_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-14.44)/(mean_distance*14.44)/1.889726)+(source_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-12.8)/(mean_distance*12.8)/1.889726)
        Lambda=0.8+alpha_reorg*(Lambda*27.2114-0.8)
        deltaG=source_node.E-target_node.E
        VBB=2.7/(source_node.AtomQuant*source_node.AtomQuant)**0.5 * np.exp(-0.72*shortest_distance)
        VDD=10**(1.73-0.42*(shortest_distance+2*1.87))
        VIF=(VDD*VBB)**0.5        

    #Not included in the study, but for broader use. Compute Marcus parameters between the electron donor node (heme) and a bridge node (amino acid) for the forward *ET* step.  
    if type_da=='electron' and (source_node.nodeNum in donor_list and source_node.Residuetype in ['HEM']) and (target_node.nodeNum in bridge_list and target_node.Residuetype in ['CYS','MET','TRP','TYR']):
        Lambda=(source_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-17.45)/(mean_distance*17.45)/1.889726)+(target_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-12.8)/(mean_distance*12.8)/1.889726)
        Lambda=0.8+alpha_reorg*(Lambda*27.2114-0.8)
        deltaG=source_node.E-target_node.E
        VBB=2.7/(target_node.AtomQuant*target_node.AtomQuant)**0.5 * np.exp(-0.72*shortest_distance)
        VDD=np.mean(np.array([8.0])*np.exp(-1*np.array([2.6])*((shortest_distance + 2*3.63)-16.5)/2))*0.00012398
        VIF=(VDD*VBB)**0.5

    #Not included in the study, but for broader use. Compute Marcus parameters between the electron donor node (heme) and a bridge node (amino acid) for the backward *ET* step.      
    if type_da=='electron' and (source_node.nodeNum in bridge_list and source_node.Residuetype in ['CYS','MET','TRP','TYR']) and (target_node.nodeNum in donor_list and target_node.Residuetype in ['HEM']):
        Lambda=(target_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-17.45)/(mean_distance*17.45)/1.889726)+(source_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-12.8)/(mean_distance*12.8)/1.889726)
        Lambda=0.8+alpha_reorg*(Lambda*27.2114-0.8)
        deltaG=source_node.E-target_node.E
        VBB=2.7/(source_node.AtomQuant*source_node.AtomQuant)**0.5 * np.exp(-0.72*shortest_distance)
        VDD=np.mean(np.array([8.0])*np.exp(-1*np.array([2.6])*((shortest_distance + 2*3.63)-16.5)/2))*0.00012398
        VIF=(VDD*VBB)**0.5   

##############Bridge-Acceptor##############

    #Compute Marcus parameters between a bridge node (amino acid) and the terminal hole acceptor node (amino acid) for the forward HT step. In our kinetic model, we assume that there is no backward HT step from the terminal hole acceptor.
    if type_da=='hole' and (source_node.nodeNum in bridge_list and source_node.Residuetype in ['CYS','MET','TRP','TYR']) and (target_node.nodeNum in acceptor_list and target_node.Residuetype in ['CYS','MET','TRP','TYR']):
        Lambda=(source_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-12.8)/(mean_distance*12.8)/1.889726)+(target_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-12.8)/(mean_distance*12.8)/1.889726)
        Lambda=0.8+alpha_reorg*(Lambda*27.2114-0.8)
        deltaG=target_node.E-source_node.E
        VIF=2.7/(source_node.AtomQuant*target_node.AtomQuant)**0.5 * np.exp(-0.72*shortest_distance)
    
    #Not included in the study, but for broader use. Compute Marcus parameters between a bridge node (amino acid) and the terminal hole acceptor node (amino acid) for the forward *ET* step.
    if type_da=='electron' and (source_node.nodeNum in bridge_list and source_node.Residuetype in ['CYS','MET','TRP','TYR']) and (target_node.nodeNum in acceptor_list and target_node.Residuetype in ['CYS','MET','TRP','TYR']):
        Lambda=(source_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-12.8)/(mean_distance*12.8)/1.889726)+(target_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-12.8)/(mean_distance*12.8)/1.889726)
        Lambda=0.8+alpha_reorg*(Lambda*27.2114-0.8)
        deltaG=source_node.E-target_node.E
        VIF=2.7/(source_node.AtomQuant*target_node.AtomQuant)**0.5 * np.exp(-0.72*shortest_distance)

    #For BSS: Compute Marcus parameters between a bridge node (amino acid) and the terminal hole acceptor node ([4Fe4S]) for the forward HT step.
    if type_da=='hole' and (source_node.nodeNum in bridge_list and source_node.Residuetype in ['CYS','MET','TRP','TYR']) and (target_node.nodeNum in acceptor_list and target_node.Residuetype in ['SF4']):
        Lambda=(target_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-14.44)/(mean_distance*14.44)/1.889726)+(source_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-12.8)/(mean_distance*12.8)/1.889726)
        Lambda=0.8+alpha_reorg*(Lambda*27.2114-0.8)
        deltaG=target_node.E-source_node.E
        VBB=2.7/(source_node.AtomQuant*source_node.AtomQuant)**0.5 * np.exp(-0.72*shortest_distance)
        VAA=10**(1.73-0.42*(shortest_distance+2*1.87))
        VIF=(VBB*VAA)**0.5
    
    #For primase-DNA: Compute Marcus parameters between a bridge node (amino acid) and the terminal *electron* acceptor node (G/A) for the forward *ET* step.     
    if type_da=='electron' and (source_node.nodeNum in bridge_list and source_node.Residuetype in ['CYS','MET','TRP','TYR']) and (target_node.nodeNum in acceptor_list and target_node.Residuetype in ['GTP','G','DG','A','DA']):
        Lambda=(source_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-12.8)/(mean_distance*12.8)/1.889726)+(target_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(1/target_node.r-1/mean_distance)*(1/1.889726))
        Lambda=0.8+alpha_reorg*(Lambda*27.2114-0.8)
        deltaG=source_node.E-target_node.E
        VIF=2.7/(source_node.AtomQuant*target_node.AtomQuant)**0.5 * np.exp(-0.72*shortest_distance)

    #Compute Marcus parameters between the hole bridge node (amino acid) and the terminal hole acceptor node (heme) for the forward HT step.
    if type_da=='hole' and (source_node.nodeNum in bridge_list and source_node.Residuetype in ['CYS','MET','TRP','TYR']) and (target_node.nodeNum in acceptor_list and target_node.Residuetype in ['HEM']):
        Lambda=(target_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-17.45)/(mean_distance*17.45)/1.889726)+(source_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-12.8)/(mean_distance*12.8)/1.889726)
        Lambda=0.8+alpha_reorg*(Lambda*27.2114-0.8)
        deltaG=target_node.E-source_node.E
        VBB=2.7/(source_node.AtomQuant*source_node.AtomQuant)**0.5 * np.exp(-0.72*shortest_distance)
        VAA=np.mean(np.array([8.0])*np.exp(-1*np.array([2.6])*((shortest_distance + 2*3.63)-16.5)/2))*0.00012398
        VIF=(VBB*VAA)**0.5

     #Compute Marcus parameters between the hole bridge node (heme) and a terminal acceptor node (heme) for the forward HT step.  
    if type_da=='hole' and (source_node.nodeNum in bridge_list and source_node.Residuetype in ['HEM']) and (target_node.nodeNum in acceptor_list and target_node.Residuetype in ['HEM']):
        Lambda=(source_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-17.45)/(mean_distance*17.45)/1.889726)+(target_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-17.45)/(mean_distance*17.45)/1.889726)
        Lambda=0.8+alpha_reorg*(Lambda*27.2114-0.8)
        deltaG=target_node.E-source_node.E
        VBB=np.mean(np.array([8.0])*np.exp(-1*np.array([2.6])*((shortest_distance + 2*3.63)-16.5)/2))*0.00012398
        VAA=np.mean(np.array([8.0])*np.exp(-1*np.array([2.6])*((shortest_distance + 2*3.63)-16.5)/2))*0.00012398
        VIF=(VAA*VBB)**0.5

##############Donor-Acceptor##############

    #Compute Marcus parameters between the hole donor node (heme) and the terminal hole acceptor node (amino acid) for the forward HT step.
    if type_da=='hole' and (source_node.nodeNum in donor_list and source_node.Residuetype in ['HEM']) and (target_node.nodeNum in acceptor_list and target_node.Residuetype in ['CYS','MET','TRP','TYR']):
        Lambda=(source_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-17.45)/(mean_distance*17.45)/1.889726)+(target_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-12.8)/(mean_distance*12.8)/1.889726)
        Lambda=0.8+alpha_reorg*(Lambda*27.2114-0.8)
        deltaG=target_node.E-source_node.E
        VAA=2.7/(target_node.AtomQuant*target_node.AtomQuant)**0.5 * np.exp(-0.72*shortest_distance)
        VDD=np.mean(np.array([8.0])*np.exp(-1*np.array([2.6])*((shortest_distance + 2*3.63)-16.5)/2))*0.00012398
        VIF=(VDD*VAA)**0.5
    
    #Not included in the study, but for broader use. Compute Marcus parameters between the electron donor node (heme) and the terminal electron acceptor node (amino acid) for the forward *ET* step.
    if type_da=='electron' and (source_node.nodeNum in donor_list and source_node.Residuetype in ['HEM']) and (target_node.nodeNum in acceptor_list and target_node.Residuetype in ['CYS','MET','TRP','TYR']):
        Lambda=(source_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-17.45)/(mean_distance*17.45)/1.889726)+(target_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-12.8)/(mean_distance*12.8)/1.889726)
        Lambda=0.8+alpha_reorg*(Lambda*27.2114-0.8)
        deltaG=source_node.E-target_node.E
        VAA=2.7/(target_node.AtomQuant*target_node.AtomQuant)**0.5 * np.exp(-0.72*shortest_distance)
        VDD=np.mean(np.array([8.0])*np.exp(-1*np.array([2.6])*((shortest_distance + 2*3.63)-16.5)/2))*0.00012398
        VIF=(VDD*VAA)**0.5

    #Compute Marcus parameters between the hole donor node (heme) and an acceptor node (heme) for the forward HT step.  
    if type_da=='hole' and (source_node.nodeNum in donor_list and source_node.Residuetype in ['HEM']) and (target_node.nodeNum in acceptor_list and target_node.Residuetype in ['HEM']):
        Lambda=(source_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-17.45)/(mean_distance*17.45)/1.889726)+(target_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-17.45)/(mean_distance*17.45)/1.889726)
        Lambda=0.8+alpha_reorg*(Lambda*27.2114-0.8)
        deltaG=target_node.E-source_node.E
        VAA=np.mean(np.array([8.0])*np.exp(-1*np.array([2.6])*((shortest_distance + 2*3.63)-16.5)/2))*0.00012398
        VDD=np.mean(np.array([8.0])*np.exp(-1*np.array([2.6])*((shortest_distance + 2*3.63)-16.5)/2))*0.00012398
        VIF=(VDD*VAA)**0.5

    #For BSS: Compute Marcus parameters between the hole donor node (Gly) and the terminal hole acceptor node ([4Fe4S]) for the forward HT step.
    if type_da=='hole' and (source_node.nodeNum in donor_list and source_node.Residuetype in ['GLY']) and (target_node.nodeNum in acceptor_list and target_node.Residuetype in ['SF4']):
        Lambda=(source_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-12.8)/(mean_distance*12.8)/1.889726)+(target_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-14.44)/(mean_distance*14.44)/1.889726)
        Lambda=0.8+alpha_reorg*(Lambda*27.2114-0.8)
        deltaG=target_node.E-source_node.E
        VDD=2.7/(source_node.AtomQuant*source_node.AtomQuant)**0.5 * np.exp(-0.72*shortest_distance)
        VAA=10**(1.73-0.42*(shortest_distance+2*1.87))
        VIF=(VDD*VAA)**0.5
        
    #For primase-DNA: Compute Marcus parameters between the *electron* donor node ([4Fe4S]) and the terminal *electron* acceptor node (G/A) for the forward *ET* step.     
    if type_da=='electron' and (source_node.nodeNum in donor_list and source_node.Residuetype in ['SF4']) and (target_node.nodeNum in acceptor_list and target_node.Residuetype in ['GTP','G','DG','A','DA']):
        Lambda=(source_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-14.44)/(mean_distance*14.44)/1.889726)+(target_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(1/target_node.r-1/mean_distance)*(1/1.889726))
        Lambda=0.8+alpha_reorg*(Lambda*27.2114-0.8)
        deltaG=source_node.E-target_node.E
        VAA=2.7/(target_node.AtomQuant*target_node.AtomQuant)**0.5 * np.exp(-0.72*shortest_distance)
        VDD=10**(1.73-0.42*(shortest_distance+2*1.87))
        VIF=(VDD*VAA)**0.5
        
##############Bridge-Bridge##############

    #Compute Marcus parameters between two bridge nodes (amino acid).  
    if type_da=='hole' and (source_node.nodeNum in bridge_list and source_node.Residuetype in ['CYS','MET','TRP','TYR']) and (target_node.nodeNum in bridge_list and target_node.Residuetype in ['CYS','MET','TRP','TYR']):
        Lambda=(source_node.Lamb+target_node.Lamb)/27.2114+(1/2.2-1/4)*(mean_distance-12.8)/(mean_distance*12.8)/1.889726
        Lambda=0.8+alpha_reorg*(Lambda*27.2114-0.8)
        deltaG=target_node.E-source_node.E
        VIF=2.7/(source_node.AtomQuant*target_node.AtomQuant)**0.5 * np.exp(-0.72*shortest_distance)
        
    if type_da=='electron' and (source_node.nodeNum in bridge_list and source_node.Residuetype in ['CYS','MET','TRP','TYR']) and (target_node.nodeNum in bridge_list and target_node.Residuetype in ['CYS','MET','TRP','TYR']):
        Lambda=(source_node.Lamb+target_node.Lamb)/27.2114+(1/2.2-1/4)*(mean_distance-12.8)/(mean_distance*12.8)/1.889726
        Lambda=0.8+alpha_reorg*(Lambda*27.2114-0.8)
        deltaG=source_node.E-target_node.E
        VIF=2.7/(source_node.AtomQuant*target_node.AtomQuant)**0.5 * np.exp(-0.72*shortest_distance)

    #Compute Marcus parameters between two bridge nodes (heme).  
    if type_da=='hole' and (source_node.nodeNum in bridge_list and source_node.Residuetype in ['HEM']) and (target_node.nodeNum in bridge_list and target_node.Residuetype in ['HEM']):
        Lambda=(source_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-17.45)/(mean_distance*17.45)/1.889726)+(target_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-17.45)/(mean_distance*17.45)/1.889726)
        Lambda=0.8+alpha_reorg*(Lambda*27.2114-0.8)
        deltaG=target_node.E-source_node.E
        VBB1=np.mean(np.array([8.0])*np.exp(-1*np.array([2.6])*((shortest_distance + 2*3.63)-16.5)/2))*0.00012398
        VBB2=np.mean(np.array([8.0])*np.exp(-1*np.array([2.6])*((shortest_distance + 2*3.63)-16.5)/2))*0.00012398
        VIF=(VBB1*VBB2)**0.5

    #Compute Marcus parameters between a bridge node (heme) and bridge node (amino acid) for the forward HT step.  
    if type_da=='hole' and (source_node.nodeNum in bridge_list and source_node.Residuetype in ['HEM']) and (target_node.nodeNum in bridge_list and target_node.Residuetype in ['CYS','MET','TRP','TYR']):
        Lambda=(source_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-17.45)/(mean_distance*17.45)/1.889726)+(target_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-12.8)/(mean_distance*12.8)/1.889726)
        Lambda=0.8+alpha_reorg*(Lambda*27.2114-0.8)
        deltaG=target_node.E-source_node.E
        VBB1=2.7/(target_node.AtomQuant*target_node.AtomQuant)**0.5 * np.exp(-0.72*shortest_distance)
        VBB2=np.mean(np.array([8.0])*np.exp(-1*np.array([2.6])*((shortest_distance + 2*3.63)-16.5)/2))*0.00012398
        VIF=(VBB1*VBB2)**0.5

    #Compute Marcus parameters between a bridge donor node (heme) and a bridge node (amino acid) for the backward HT step.      
    if type_da=='hole' and (source_node.nodeNum in bridge_list and source_node.Residuetype in ['CYS','MET','TRP','TYR']) and (target_node.nodeNum in bridge_list and target_node.Residuetype in ['HEM']):
        Lambda=(target_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-17.45)/(mean_distance*17.45)/1.889726)+(source_node.Lamb/27.2114+0.5*(1/2.2-1/4)*(mean_distance-12.8)/(mean_distance*12.8)/1.889726)
        Lambda=0.8+alpha_reorg*(Lambda*27.2114-0.8)
        deltaG=target_node.E-source_node.E
        VBB1=2.7/(source_node.AtomQuant*source_node.AtomQuant)**0.5 * np.exp(-0.72*shortest_distance)
        VBB2=np.mean(np.array([8.0])*np.exp(-1*np.array([2.6])*((shortest_distance + 2*3.63)-16.5)/2))*0.00012398
        VIF=(VBB1*VBB2)**0.5

#Nonadiabatic expression for the HT rate

    smallK= ((2*3.141592654)**2/4.1357e-15)*(VIF**2)*(4*3.141592654*Lambda*298.15*8.6173303e-5)**-0.5*np.exp(-1*(deltaG+Lambda)**2/(4*Lambda*8.6173303e-5*298.15))     
    return smallK

##############Defining edges in directed graph##############
donor_bridge_combine=[(x,y) for x in donor_list for y in bridge_list]
import itertools
bridge_combine=list(itertools.combinations(bridge_list, 2))
bridge_acceptor_combine=[(x,y) for x in bridge_list for y in acceptor_list]
donor_acceptor_combine=[(x,y) for x in donor_list for y in acceptor_list]                  

#Add edges between nodes.
for e in donor_bridge_combine:
    G.add_edge(e[0],e[1],weight=smallK(e[0],e[1]))
    G.add_edge(e[1],e[0],weight=smallK(e[1],e[0]))

for e in bridge_combine:
    G.add_edge(e[0],e[1],weight=smallK(e[0],e[1]))
    G.add_edge(e[1],e[0],weight=smallK(e[1],e[0]))

#One directed edge is considered for bridge-acceptor or donor-acceptor as the backward HT step is assumed to be negligible.
for e in bridge_acceptor_combine:
    G.add_edge(e[0],e[1],weight=smallK(e[0],e[1]))

for e in donor_acceptor_combine:
    G.add_edge(e[0],e[1],weight=smallK(e[0],e[1]))

def get_k_list(path):
    edge_list=[(path[i],path[i+1]) for i in range(len(path)-1)]
    k_list=[G[e[0]][e[1]]['weight'] for e in edge_list]
    return k_list

##############Define the approximate expression for the mean residence time (Equation 7)##############
from numpy import prod 
def bigK(k_list):
    return 1.0/sum([1/e for e in k_list])

def calculate_all_K(all_pathways):
    all_k_lists=[get_k_list(e) for e in all_pathways]
    all_big_K=[bigK(e) for e in all_k_lists]
    result={}
    for i in range(len(all_pathways)):
        result[str(all_pathways[i])]=all_big_K[i]
    return result

all_big_K={}
for e in acceptor_list:
    paths = nx.all_simple_paths(G,source=donor_num,target=e,cutoff=cutoff_num)
    all_possible_paths=list(paths)
    temp=calculate_all_K(all_possible_paths)
    all_big_K= {**all_big_K,**temp}
print(' ')
print("Top" + " {} ".format(total_paths)  + "pathways ranked according to the approximate mean residence time:")
#Top X pathways ranked according to the approximate expression for the mean residence time.
for i in range(total_paths):
    print(max(all_big_K, key=all_big_K.get),1/all_big_K[max(all_big_K, key=all_big_K.get)])
    s=max(all_big_K, key=all_big_K.get)
    all_big_K.pop(s)

##############Define the exact expression for the mean residence time (Equation 6)##############
def bigK_hard(path): 
    tau=0
    N=len(path)-1 
    for n in range(0,N):
        q=0
        
        for j in range(0,N-n-1):
            p=1            
            for i in range(n+1,N-j):
                s=G[path[i]][path[i-1]]['weight']/G[path[i]][path[i+1]]['weight']
                p=p*s
            q=q+p
        tau=tau+1/G[path[n]][path[n+1]]['weight']*(q+1)
    return 1/tau

def calculate_all_K_hard(all_pathways): 
    all_big_K=[bigK_hard(e) for e in all_pathways]
    result={}
    for i in range(len(all_pathways)):
        result[str(all_pathways[i])]=all_big_K[i]
    return result

all_big_K_hard={}
for e in acceptor_list:
    paths = nx.all_simple_paths(G, source=donor_num, target=e,cutoff=cutoff_num) 
    all_possible_paths=list(paths)
    temp=calculate_all_K_hard(all_possible_paths)
    all_big_K_hard= {**all_big_K_hard,**temp}
print(' ')
print("Top" + " {} ".format(total_paths)  + "pathways ranked according to the exact mean residence time:")
#Top X pathways ranked according to the exact expression for the mean residence time.
for i in range(total_paths):
    print(max(all_big_K_hard, key=all_big_K_hard.get),1/all_big_K_hard[max(all_big_K_hard, key=all_big_K_hard.get)])
    s=max(all_big_K_hard, key=all_big_K_hard.get)
    all_big_K_hard.pop(s)
print(' ')
print('End of calculation')
print(' ')

##############Section for allowing user to obtain the exact mean residence time of a specific pathway##############
path_calculate=input('Would you like to find the exact mean residence time of a specific pathway? (Yes/No): ')
print(' ')
import time 
if path_calculate=='Yes': 
        pathsa=input('Please specify the pathway with their residue numbers (for example, 999, 156, 115, 112, 305): ')
        time.sleep(1)
        pathsb=[[int(s) for s in pathsa.split(',')]]
        all_possible_paths=list(pathsb)
        all_big_K_hard=calculate_all_K_hard(all_possible_paths)
        print(' ')
        print('Exact mean residence time: ',1/all_big_K_hard[max(all_big_K_hard, key=all_big_K_hard.get)])
print(' ')
time.sleep(3)
print('*****************************************')
print('If you have any questions or suggestions about EHPath, please contact Ruijie Teo (rt131@duke.edu) or Dr. Ruobing Wang (ruobing.wang@operasolutions.com).')
print('Citation: Teo, R. D.; Wang, R.; Smithwick, E.; Migliore, A.; Therien, M. J.; Beratan, D. N. Proc. Natl. Acad. Sci. U.S.A., 2019, 116, 15811-15816.')
print('*****************************************')
print(' ')
print('Copyright (C) 2019  Teo, R. D.; Wang, R.; Smithwick, E.; Migliore, A.; Therien, M. J.; Beratan, D. N.')
print(' ')
print('This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.')
print(' ')
print('This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.')
print(' ')
print('You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/gpl.html.')
print(' ')
print('Press Enter to exit EHPath.')
input()
