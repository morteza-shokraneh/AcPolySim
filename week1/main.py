#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  7 01:20:19 2022

@author: morteza
"""

import numpy as np
import matplotlib.pyplot as plt
import sys 

###############################################################################
class Polymer:
    def __init__(self, network, L):
        self.L = L
        self.network = network
        self.confs = np.zeros((2, L+1), dtype=np.int64)

        while True:
            x0 = np.random.randint(0, len(network)); y0 = np.random.randint(0, len(network))
            if network[y0, x0] == 0:
                self.confs[0, 0] = y0
                self.confs[1, 0] = x0
                self._fill_the_network(0, is_grow=True, is_head=False)
                break
            else:
                continue

    def _fill_the_network(self, *args, is_grow: bool, is_head: bool):
        if is_grow == True:
            self.network[self.confs[0, args], self.confs[1, args]] = 1
        elif is_grow == False and is_head == True:
            self.network[self.confs[0, self.L+1], self.confs[1, self.L+1]] = 1
            self.network[self.confs[0, 0], self.confs[1, 0]] = 0

        elif is_grow == False and is_head == False:
            self.network[self.confs[0, 0], self.confs[1, 0]] = 1
            self.network[self.confs[0, self.L+1], self.confs[1, self.L+1]] = 0

    def _implement_boundary_condition(self, index, move):
        return (self.confs[0:,index]+move)%(len(self.network))   

    def _find_open_site(self, pos):
        move_array = np.array([[1, 0], [-1, 0], [0, 1], [0, -1]], dtype=np.int64)
        success = False
        while success == False:
            if move_array.shape[0] == 0:
                return self.confs[0:, pos]
            index = np.random.randint(0, move_array.shape[0])
            move = move_array[index, 0:]
            possible_destination = self._implement_boundary_condition(pos, move)
            if self.network[possible_destination[0], possible_destination[1]] == 0:
                return possible_destination 
            else:
                move_array = np.delete(move_array, index, axis=0)

    def _generate_single_chain(self):
        for step in range(1, self.L+1):
            destination = self._find_open_site(step-1)
            if np.array_equal(self.confs[0:, step-1], destination):
                self.confs = np.vstack((np.trim_zeros(self.confs[0, 0:]), np.trim_zeros(self.confs[1, 0:])))
                self.L = self.confs.shape[1] - 1
                print("the road is blocked!\n")
                print("the length of generated chain is",len(np.trim_zeros(self.confs[1, 0:])))
                break
            else:
                self.confs[0:, step] = destination
                self._fill_the_network(step, is_grow=True, is_head=True)

        return self.confs
    
    def _grow(self, is_head:bool):
        if is_head==True:
            pos=self.L 
            destination = self._find_open_site(pos) 
            if np.array_equal(self.confs[0:,pos], destination):
                return self._grow(is_head=False)
                #return print("the update cannot proceed") 
            else:
                self.confs = np.concatenate((self.confs,destination.reshape(2,1)),axis=1)
                self.L+=1; pos=self.L 
                self._fill_the_network(pos, is_grow=True,is_head=True) 
        else:
            destination = self._find_open_site(0)
            if np.array_equal(self.confs[0:,0], destination):
                return self._grow(is_head=True)
                #return print("the update cannot proceed")
            else:
                self.confs = np.concatenate((destination.reshape(2,1),self.confs),axis=1)
                self.L+=1; pos=0
                self._fill_the_network(pos, is_grow=True, is_head=False)
       
                
    def _move(self, is_head:bool):
        if is_head==True:
            pos=self.L
            destination = self._find_open_site(pos)
            if np.array_equal(self.confs[0:,pos], destination):
                #return print("the update cannot proceed")
                return self._move(is_head=False)
                #return sys.exit("Error")
            else:
                self.confs = np.concatenate((self.confs,destination.reshape(2,1)),axis=1)
                self._fill_the_network(is_grow=False, is_head=True)
                self.confs = np.delete(self.confs,0,axis=1)
        else:
            pos=0
            destination = self._find_open_site(pos) 
            if np.array_equal(self.confs[0:,0], destination):
                return self._move(is_head=True)
                #print("the update cannot proceed")
            else:
                self.confs = np.concatenate((destination.reshape(2,1),self.confs),axis=1)
                self._fill_the_network(is_grow=False, is_head=False)
                self.confs = np.delete(self.confs, self.L+1, axis=1)
            
                

        
    def update(self,gr_pr,ht_pr):
        if np.random.rand()<=ht_pr:
            self._grow(is_head=True) if np.random.rand()<=gr_pr else self._move(is_head=True)
        else:
            self._grow(is_head=False) if np.random.rand()<=gr_pr else self._move(is_head=False)
            


    def visualize(self, time):
        confs = np.vstack((np.trim_zeros(self.confs[0, 0:]), np.trim_zeros(self.confs[1, 0:])))
        where = 0
        for site in range(confs.shape[1]):
            if confs[0, site] == len(self.network) and confs[0, site+1] == 0:
                plt.plot(confs[1, where:site+1],confs[0, where:site+1], color='blue')
                where = site + 1
            elif confs[0, site] == 0 and confs[0, site+1] == len(self.network):
                plt.plot(confs[1, where:site+1],confs[0, where:site+1], color='blue') 
                where = site + 1
            elif confs[1, site] == len(self.network) and confs[1, site+1] == 0:
                plt.plot(confs[1, where:site+1],confs[0, where:site+1], color='blue')
                where = site + 1
            elif confs[1, site] == 0 and confs[1, site+1] == len(self.network):
                plt.plot(confs[1,where:site+1],confs[0,where:site+1],color='blue')
                where = site + 1
            else:
                plt.plot(confs[1,where:site+1],confs[0,where:site+1],color='blue')
        print("the length of current chain is", confs.shape[1])       
        plt.legend()
        if time<10:
            plt.savefig("test1-40-70/"+str(0)+str(0)+str(time)+'.png') 
        elif time>=10 and time<100:
            plt.savefig("test1-40-70/"+str(0)+str(time)+'.png')
        else:
            plt.savefig("test1-40-70/"+str(time)+'.png')  
        plt.show() 
###############################################################################


network = np.zeros((60,60),dtype=np.int64) 
L = 40
test = Polymer(network, L) 
test._generate_single_chain()

 
for i in range(300):
    test.update(0.4, 0.7)
    print(test.confs)
    test.visualize(i)
