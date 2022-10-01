#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 04:09:43 2022

@author: morteza shokraneh
"""
import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

################## Parameters #################################################
L = 10
N = 20
size = 10
ht_pr = 0.7
gr_pr = 0.6
###############################################################################



################ Three D Polymer Class ########################################
class Polymer:
    def __init__(self, network, L):
        self.original_network = np.copy(network) 
        self.network = np.copy(network)  
        self.L = L 
        self._generate_single_chain() 
            
                        
    def _implement_boundary_condition(self, index, move):
        return (self.confs[0:,index]+move)%len(self.network)
    
    
    def _fill_the_network(self, *args, is_grow:bool, is_head:bool, is_generate:bool):
        if is_generate==True:
            self.network[self.confs[0,args],self.confs[1,args],self.confs[2,args]] = 1
        elif is_grow==True and is_head==True and is_generate==False:
            self.network[self.confs[0,self.L],self.confs[1,self.L],self.confs[2,self.L]] = 1
        elif is_grow==True and is_head==False and is_generate==False:
            self.network[self.confs[0,0],self.confs[1,0],self.confs[2,0]] = 1
        elif is_grow==False and is_head==True and is_generate==False:
            self.network[self.confs[0,self.L+1],self.confs[1,self.L+1],self.confs[2,self.L+1]] = 1
            self.network[self.confs[0,0],self.confs[1,0],self.confs[2,0]] = 0
        elif is_grow==False and is_head==False and is_generate==False:
            self.network[self.confs[0,0],self.confs[1,0],self.confs[2,0]] = 1
            self.network[self.confs[0,self.L+1],self.confs[1,self.L+1],self.confs[2,self.L+1]] = 0
            
            
    def calculate_acceptance_rate(self):
        pos = 0
        move_array = np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]], dtype=np.int64)
        self.acceptance_rate = 0
        for index in range(move_array.shape[0]):
            possible_destination = self._implement_boundary_condition(pos, move_array[index, 0:])
            if self.network[possible_destination[0], possible_destination[1], possible_destination[2]]==0:
                self.acceptance_rate += (1/move_array.shape[0])
        pos = self.L 
        move_array = np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]], dtype=np.int64)
        #self.acceptance_rate = 0
        for index in range(move_array.shape[0]):
            possible_destination = self._implement_boundary_condition(pos, move_array[index, 0:])
            if self.network[possible_destination[0], possible_destination[1], possible_destination[2]]==0:
                self.acceptance_rate += (1/move_array.shape[0])
        self.acceptance_rate = self.acceptance_rate / 2.0
        return self.acceptance_rate
    


                    
    def _find_open_site(self, pos):
        move_array = np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]], dtype=np.int64)
        success = False
        while success==False:
            if move_array.shape[0]==0:
                return self.confs[0:, pos]
            index = np.random.randint(0, move_array.shape[0])
            move = move_array[index, 0:]
            possible_destination = self._implement_boundary_condition(pos, move)
            if self.network[possible_destination[0], possible_destination[1],possible_destination[2]]==0:
                return possible_destination
            else:
                move_array = np.delete(move_array, index, axis=0)
            
                        
    def _generate_single_chain(self):
        self.confs = np.zeros((3,self.L+1), dtype=np.int64)
        while True:
            x0 = np.random.randint(0, len(self.network))
            y0 = np.random.randint(0, len(self.network))
            z0 = np.random.randint(0, len(self.network))
            if self.network[x0,y0,z0] == 0:
                self.confs[0,0] = x0; self.confs[1,0]=y0; self.confs[2,0]=z0
                self._fill_the_network(0, is_grow=True, is_head=True, is_generate=True)
                break
            else:
                continue
            
        for step in range(self.L):
            destination = self._find_open_site(step)
            if np.array_equal(self.confs[0:,step], destination)==True:
                print("attempt failed, trying again...")
                self.network = np.copy(self.original_network)
                return self._generate_single_chain()
            else:
                self.confs[0:,step+1]=destination
                self._fill_the_network(step+1, is_grow=True, is_head=True, is_generate=True)
 
                
    def _grow(self, is_head:bool):
        if is_head==True:
            pos = self.L
            #self.calculate_acceptance_rate(pos)
            destination = self._find_open_site(pos)
            if np.array_equal(self.confs[0:,pos], destination)==True:
                return None 
                #return self._grow(is_head=False)
            else:
                self.confs = np.concatenate((self.confs,destination.reshape(3,1)),axis=1)
                self.L+=1; pos=self.L
                self._fill_the_network(is_grow=True, is_head=True, is_generate=False)
        else:
            pos = 0
            #self.calculate_acceptance_rate(pos)
            destination = self._find_open_site(pos)
            if np.array_equal(self.confs[0:,0], destination)==True:
                return None 
                #return self._grow(is_head=True)
            else:
                self.confs = np.concatenate((destination.reshape(3,1),self.confs),axis=1)
                self.L+=1; pos=0;
                self._fill_the_network(is_grow=True, is_head=False, is_generate=False)
                
                
    def _move(self, is_head:bool):
        if is_head==True:
            pos=self.L
            #self.calculate_acceptance_rate(pos)
            destination = self._find_open_site(pos)
            if np.array_equal(self.confs[0:,pos], destination)==True:
                return None 
                #return self._move(is_head=False)
            else:
                self.confs = np.concatenate((self.confs,destination.reshape(3,1)),axis=1)
                self._fill_the_network(is_grow=False, is_head=True, is_generate=False)
                self.confs = np.delete(self.confs,0,axis=1)
        else:
            pos = 0
            #self.calculate_acceptance_rate(pos)
            destination = self._find_open_site(pos)
            if np.array_equal(self.confs[0:,pos], destination)==True:
                return None 
                #self._move(is_head=True)
            else:
                self.confs = np.concatenate((destination.reshape(3,1),self.confs),axis=1)
                self._fill_the_network(is_grow=False, is_head=False, is_generate=False)
                self.confs = np.delete(self.confs, self.L+1, axis=1)
                
                
    def update(self, ht_pr, gr_pr):
        if np.random.rand() <= ht_pr:
            self._grow(is_head=True) if np.random.rand()<=gr_pr else self._move(is_head=False)
        else:
            self._grow(is_head=False) if np.random.rand()<=gr_pr else self._move(is_head=False)
            
            
    def update_network(self, network):
        self.network = np.copy(network)


    def visualize(self, color, ax):
        where = 0
        for point in range(self.L):
            if self.confs[2,point]==0 and self.confs[2,point+1]==9:
                ax.plot3D(self.confs[0,where:(point+1)],self.confs[1,where:(point+1)],self.confs[2,where:(point+1)],color=color)
                where=point+1
            elif self.confs[2,point]==9 and self.confs[2,point+1]==0:
                ax.plot3D(self.confs[0,where:(point+1)],self.confs[1,where:(point+1)],self.confs[2,where:(point+1)],color=color)
                where=point+1
            elif self.confs[1,point]==0 and self.confs[1,point+1]==9:
                ax.plot3D(self.confs[0,where:(point+1)],self.confs[1,where:(point+1)],self.confs[2,where:(point+1)],color=color)
                where=point+1
            elif self.confs[1,point]==9 and self.confs[1,point+1]==0:
                ax.plot3D(self.confs[0,where:(point+1)],self.confs[1,where:(point+1)],self.confs[2,where:(point+1)],color=color)
                where=point+1
            elif self.confs[0,point]==0 and self.confs[0,point+1]==9:
                ax.plot3D(self.confs[0,where:(point+1)],self.confs[1,where:(point+1)],self.confs[2,where:(point+1)],color=color)
                where=point+1
            elif self.confs[0,point]==9 and self.confs[0,point+1]==0:
                ax.plot3D(self.confs[0,where:(point+1)],self.confs[1,where:(point+1)],self.confs[2,where:(point+1)],color=color)
                where=point+1
        ax.plot3D(self.confs[0,where:(self.L+1)],self.confs[1,where:(self.L+1)],self.confs[2,where:(self.L+1)],color=color) 
######################  END  ##################################################       
        
                
####################  Three D Melt Class  #####################################
class Melt:
    def __init__(self, N, L, size):
        self.network = np.zeros((size,size,size),dtype=np.int64)
        self.L = L
        self.N = N
        self.polymers = [] 
        self.acceptance_rate = 0
        
        for poly in range(N):
            local_polymer = Polymer(self.network, L)
            self.polymers.append(local_polymer)
            self.network = local_polymer.network
                       
        self.colors=[]
        for poly in range(N):
            self.colors.append('#%06X' % np.random.randint(0, 0xFFFFFF))
            
    def _get_acceptance_rate(self):
        for poly in range(self.N):
            self.acceptance_rate += self.polymers[poly].calculate_acceptance_rate()
        self.acceptance_rate = self.acceptance_rate / self.N 
            
    
    def update(self, ht_pr, gr_pr):
        self._get_acceptance_rate()
        poly = np.random.randint(0,self.N)
        self.polymers[poly].update_network(self.network)
        self.polymers[poly].update(ht_pr, gr_pr)
        #self.acceptance_rate = self.polymers[poly].acceptance_rate
        self.network = self.polymers[poly].network
        

    def info(self):
        self.volume_fraction = 0
        for poly in range(self.N):
            self.volume_fraction += (self.polymers[poly].confs.shape[1]-1)
        print("current volume fraction is", self.volume_fraction)
        

    def visualize(self, time):
        ax = plt.axes(projection='3d')
        for poly in range(self.N):
            self.polymers[poly].visualize(self.colors[poly], ax)
        ax.legend()
        #plt.show()
        self.volume_fraction = 0
        for poly in range(self.N):
            self.volume_fraction += (self.polymers[poly].confs.shape[1]-1)
        print("current volume fraction is", self.volume_fraction)
        if time<10:
            plt.savefig("result/"+str(0)+str(0)+str(time)+'.png', dpi=1200) 
        elif time>=10 and time<100:
            plt.savefig("result/"+str(0)+str(time)+'.png', dpi=1200)
        else:
            plt.savefig("result/"+str(time)+'.png', dpi=1200)
######################## END  #################################################  


#vf = np.array([])
#test = Melt(N, L, size)
#for i in range(5000):
    #test.update(ht_pr, gr_pr)
    #test.info()
    #vf = np.append(vf, test.volume_fraction)
#np.save("volume_fraction.npy", vf)


vf = np.zeros(3000)
acceptance_rate = np.zeros(3000)

for run in range(500):
    test = Melt(N, L, size)
    local_vf = np.array([])
    local_acceptance_rate = np.array([])
    for i in range(3000):
        test.update(ht_pr, gr_pr)
        test.info()
        local_vf = np.append(local_vf, test.volume_fraction)
        local_acceptance_rate = np.append(local_acceptance_rate, test.acceptance_rate)
        
    vf = vf + local_vf
    acceptance_rate = acceptance_rate + local_acceptance_rate
    
    
plt.plot(vf, acceptance_rate)   
np.save("volume_fractionII.npy", vf/500)
np.save("acceptance_rate.npy", acceptance_rate/500)