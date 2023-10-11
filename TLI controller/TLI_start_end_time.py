# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#%%
start_of_bt_change = 83.45634563456346
start_of_bt_equals_zero = 125.14851485148515
end_of_bt_equals_zero = 139.40594059405942
#%%
duration_maintain = start_of_bt_equals_zero - start_of_bt_change
duration_suppress = end_of_bt_equals_zero - start_of_bt_equals_zero
#%%
maintain = 56 * (duration_maintain / (duration_maintain + duration_suppress))
suppress = 56 * (duration_suppress / (duration_maintain + duration_suppress))