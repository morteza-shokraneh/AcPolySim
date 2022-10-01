#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 10:26:56 2022

@author: morteza
"""

import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt


################## Parameters #################################################
L = 5
N = 1
L_lim = 10
size = 40
ht_pr = 0.7
gr_pr = 0.6
##############################################################################


################ Two D Polymer Class ##########################################
class TwoDPolymer:
    def __init__(self, network, L, *args, is_confs:bool):
        if is_confs==True:
            self.original_network = np.copy(network)
            self.network = np.copy(network)
            self.L = L
            self.confs = np.zeros((2,L+1))
        else:
            self.original_network = np.copy(network) 
            self.network = np.copy(network)  
            self.L = L 
            self._generate_single_chain() 
            
                        
    def _implement_boundary_condition(self, index, move):
        return (self.confs[0:,index]+move)%len(self.network)
    
    
    def _fill_the_network(self, *args, is_grow:bool, is_head:bool, is_generate:bool):
        if is_generate==True:
            self.network[self.confs[0,args],self.confs[1,args]] = 1
        elif is_grow==True and is_head==True and is_generate==False:
            self.network[self.confs[0,self.L],self.confs[1,self.L]] = 1
        elif is_grow==True and is_head==False and is_generate==False:
            self.network[self.confs[0,0],self.confs[1,0]] = 1
        elif is_grow==False and is_head==True and is_generate==False:
            self.network[self.confs[0,self.L+1],self.confs[1,self.L+1]] = 1
            self.network[self.confs[0,0],self.confs[1,0]] = 0
        elif is_grow==False and is_head==False and is_generate==False:
            self.network[self.confs[0,0],self.confs[1,0]] = 1
            self.network[self.confs[0,self.L+1],self.confs[1,self.L+1]] = 0
            
            
    def calculate_acceptance_rate(self):
        pos = 0
        move_array = np.array([[1,0],[-1,0],[0,1],[0,-1]], dtype=np.int64)
        self.acceptance_rate = 0
        for index in range(move_array.shape[0]):
            possible_destination = self._implement_boundary_condition(pos, move_array[index, 0:])
            if self.network[possible_destination[0], possible_destination[1]]==0:
                self.acceptance_rate += (1/move_array.shape[0])
        pos = self.L 
        move_array = np.array([[1,0],[-1,0],[0,1],[0,-1]], dtype=np.int64) 
        for index in range(move_array.shape[0]):
            possible_destination = self._implement_boundary_condition(pos, move_array[index, 0:])
            if self.network[possible_destination[0], possible_destination[1]]==0:
                self.acceptance_rate += (1/move_array.shape[0])
        self.acceptance_rate = self.acceptance_rate / 2.0
        return self.acceptance_rate
    


                    
    def _find_open_site(self, pos):
        move_array = np.array([[1,0],[-1,0],[0,1],[0,-1]], dtype=np.int64)
        success = False
        while success==False:
            if move_array.shape[0]==0:
                return self.confs[0:, pos]
            index = np.random.randint(0, move_array.shape[0])
            move = move_array[index, 0:]
            possible_destination = self._implement_boundary_condition(pos, move)
            if self.network[possible_destination[0], possible_destination[1]]==0:
                return possible_destination
            else:
                move_array = np.delete(move_array, index, axis=0)
            
                        
    def _generate_single_chain(self):
        self.confs = np.zeros((2,self.L+1), dtype=np.int64)
        while True:
            x0 = np.random.randint(0, len(self.network))
            y0 = np.random.randint(0, len(self.network))
            if self.network[x0,y0] == 0:
                self.confs[0,0] = x0; self.confs[1,0]=y0
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
            destination = self._find_open_site(pos)
            if np.array_equal(self.confs[0:,pos], destination)==True:
                return None 
                #return self._grow(is_head=False)
            else:
                self.confs = np.concatenate((self.confs,destination.reshape(2,1)),axis=1)
                self.L+=1; pos=self.L
                self._fill_the_network(is_grow=True, is_head=True, is_generate=False)
        else:
            pos = 0
            destination = self._find_open_site(pos)
            if np.array_equal(self.confs[0:,0], destination)==True:
                return None 
                #return self._grow(is_head=True)
            else:
                self.confs = np.concatenate((destination.reshape(2,1),self.confs),axis=1)
                self.L+=1; pos=0;
                self._fill_the_network(is_grow=True, is_head=False, is_generate=False)
            
                
                
    def _move(self, is_head:bool):
        if is_head==True:
            pos=self.L
            destination = self._find_open_site(pos)
            if np.array_equal(self.confs[0:,pos], destination)==True:
                return None 
                #return self._move(is_head=False)
            else:
                self.confs = np.concatenate((self.confs,destination.reshape(2,1)),axis=1)
                self._fill_the_network(is_grow=False, is_head=True, is_generate=False)
                self.confs = np.delete(self.confs,0,axis=1)
        else:
            pos = 0
            destination = self._find_open_site(pos)
            if np.array_equal(self.confs[0:,pos], destination)==True:
                return None 
                #self._move(is_head=True)
            else:
                self.confs = np.concatenate((destination.reshape(2,1),self.confs),axis=1)
                self._fill_the_network(is_grow=False, is_head=False, is_generate=False)
                self.confs = np.delete(self.confs, self.L+1, axis=1)
                
                
    def update(self, ht_pr, gr_pr):
        if np.random.rand() <= ht_pr:
            self._grow(is_head=True) if np.random.rand()<=gr_pr else self._move(is_head=False)
        else:
            self._grow(is_head=False) if np.random.rand()<=gr_pr else self._move(is_head=False)
    
    def update_confs(self, confs):
        self.confs = np.copy(confs)
            
            
    def update_network(self, network):
        self.network = np.copy(network)


    def visualize(self, color):
        where = 0
        for point in range(self.L):
            if self.confs[1,point]==0 and self.confs[1,point+1]==(len(self.network)-1):
                plt.plot(self.confs[0,where:(point+1)],self.confs[1,where:(point+1)],color=color)
                where=point+1
            elif self.confs[1,point]==(len(self.network)-1) and self.confs[1,point+1]==0:
                plt.plot(self.confs[0,where:(point+1)],self.confs[1,where:(point+1)],color=color)
                where=point+1
            elif self.confs[0,point]==0 and self.confs[0,point+1]==(len(self.network)-1):
                plt.plot(self.confs[0,where:(point+1)],self.confs[1,where:(point+1)],color=color)
                where=point+1
            elif self.confs[0,point]==(len(self.network)-1) and self.confs[0,point+1]==0:
                plt.plot(self.confs[0,where:(point+1)],self.confs[1,where:(point+1)],color=color) 
                where=point+1
        plt.plot(self.confs[0,where:(self.L+1)],self.confs[1,where:(self.L+1)],color=color) 
######################  END  ##################################################       
        
                
####################  Two D Melt Class  #######################################
class TwoDMelt:
    def __init__(self, N, L, size, L_lim):
        self.network = np.zeros((size,size),dtype=np.int64)
        self.L = L
        self.N = N
        self.polymers = [] 
        self.acceptance_rate = 0
        self.L_lim = L_lim
        
        for poly in range(N):
            local_polymer = TwoDPolymer(self.network, L, is_confs=False)
            self.polymers.append(local_polymer)
            self.network = local_polymer.network 
                       
        self.colors=[]
        for poly in range(1000):
            self.colors.append('#%06X' % np.random.randint(0, 0xFFFFFF))
            
    def _get_acceptance_rate(self):
        for poly in range(self.N):
            self.acceptance_rate += self.polymers[poly].calculate_acceptance_rate()
        self.acceptance_rate = self.acceptance_rate / self.N 
            
    
    def update(self, ht_pr, gr_pr):
        poly = np.random.randint(0,self.N)
        self.polymers[poly].update_network(self.network)
        self.polymers[poly].update(ht_pr, gr_pr)
        self.network = np.copy(self.polymers[poly].network)
        if self.polymers[poly].L==self.L_lim:
            self.network[self.polymers[poly].confs[0,int(self.L_lim/2)], self.polymers[poly].confs[1,int(self.L_lim/2)]] = 0
            local_polymer = TwoDPolymer(np.copy(self.network), int(self.L_lim/2 - 1),is_confs=True)
            local_polymer.update_confs(self.polymers[poly].confs[0:,:int(self.L_lim/2)])
            self.polymers.append(local_polymer)
            self.N+=1
            self.polymers[poly].confs = np.delete(self.polymers[poly].confs, np.arange(int(self.L_lim/2)+1), axis=1)
            self.polymers[poly].L = int(self.L_lim/2 - 1)
        else:
            pass
            
        
        

    def info(self):
        self.volume_fraction = 0
        for poly in range(self.N):
            self.volume_fraction += (self.polymers[poly].confs.shape[1]-1)
        print("current volume fraction is", self.volume_fraction)
        print("number of total polymer currently is", self.N)
        

    def visualize(self, time):
        for poly in range(self.N):
            self.polymers[poly].visualize(self.colors[poly])
        if time<10:
            plt.savefig("result/"+str(0)+str(0)+str(time)+'.png') 
        elif time>=10 and time<100:
            plt.savefig("result/"+str(0)+str(time)+'.png')
        else:
            plt.savefig("result/"+str(time)+'.png')
        plt.figure() 
######################## END  #################################################




test = TwoDMelt(N, L, size, L_lim)
for i in range(17000):
    test.update(ht_pr, gr_pr)
    test.info()
    if i%20==0:
        test.visualize(int(i/20))
















