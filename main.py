#-------------------------------------------------------------------------------
# Name:        SIME
# Purpose:     Scaffold-based enumeration method/software to design in-silico macrolide libraries

#
# Author:      zinph
#
# Created:     27/02/2019
# Copyright:   (c) kzphy 2019
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import time

def convert_time(second):
    '''
    Convert second to DAY, HOUR, MIN, SEC format. Parameter is second.
    '''
    day = second/86400
    hour = (day - int(day))*24
    minute = (hour - int(hour))*60
    second = round((minute - int(minute))*60,4)
    return(str(int(day)) + ' DAYS: '+ str(int(hour)) + ' HOURS: '+ str(int(minute)) + ' MINUTES: ' + str(second) + ' SECONDS')

from SIME import *

#ERY_core = 'CC[C@H]1OC(=O)[*][*sugar*][*][*sugar*][*]C[*]C(=O)[*][C@@H](O)[*]1'
ERY_core = 'CC[C@H]1OC(=O)[*][*sugar*][C@H](C)[*sugar*][*]C[*]C(=O)[*][C@@H](O)[*]1'

sample = SIME()
start_time = time.time()
sample.generate_templates_withExtendersNSugars(ERY_core)
duration = convert_time(time.time()-start_time)
print('Time Elapsed for Enumeration: ' + str(duration))
