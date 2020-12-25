# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 13:50:20 2018

@author: SUNSHINE_MACHINE
"""

from datetime import datetime

class ProgressTimer():
    def __init__(self):
        self.start_time = datetime.now()
        self.dt_list = []
        self.dt_avg = self.start_time - self.start_time
        self.approx_process_duration = None
        self.approx_process_end = None
        self.percentage = None
        
        self.n_iterations = None
        
        self._last_tick = datetime.now()
    
    def start(self,n_iterations):
        self.n_iterations = n_iterations
        self.percentage = 0
        
        self.start_time = datetime.now()
        self._last_tick = self.start_time
        
    def tick_dt(self):
        tmp_now = datetime.now()
        self.dt_list.append( tmp_now - self._last_tick )
        self._last_tick = tmp_now
        
        self._update_data()
        
    def _update_data(self):
        n = len(self.dt_list)
        
        # update average delta_t
        self.dt_avg = self.dt_avg*(float(n)-1)/n + self.dt_list[-1]/n
        
        # update approx process durability
        self.approx_process_duration = self.n_iterations*self.dt_avg
        
        # update approx time to end
        self.approx_time_to_end = self.approx_process_duration - n*self.dt_avg

        # update approx process end time
        self.approx_process_end = self.start_time + self.approx_process_duration
        
        # percentage update
        self.percentage = float(n)/self.n_iterations
        