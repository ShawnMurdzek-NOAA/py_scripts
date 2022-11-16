"""
Read a BUFR file using ecCodes

Based on this example: https://confluence.ecmwf.int/display/ECC/bufr_keys_iterator

shawn.s.murdzek@noaa.gov
Date Created: 7 October 2022
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import eccodes as ec


#---------------------------------------------------------------------------------------------------
# Decode BUFR File
#---------------------------------------------------------------------------------------------------

fname = '/mnt/lfs4/BMC/wrfruc/murdzek/sample_real_obs/obs_rap/2021072412.rap_e.t12z.prepbufr.tm00'

fptr = open(fname, 'rb')

ibufr = ec.codes_bufr_new_from_file(fptr)

# Loop through each message
n = 1
while n == 1:
#while ibufr != None:

    print()
    print('BUFR message %d' % n)
  
    try:

        # Expand message descriptors
        ec.codes_set(ibufr, 'unpack', 1)

        # Loop over each key in the message
        iterator = ec.codes_bufr_keys_iterator_new(ibufr)
        while ec.codes_bufr_keys_iterator_next(iterator):
            key = ec.codes_bufr_keys_iterator_get_name(iterator)
            print(key)  

        # Free memory
        ec.codes_release(ibufr)

        # Get next BUFR message
        ibufr = ec.codes_bufr_new_from_file(fptr)
        n = n + 1

    except ec.CodesInternalError:
    
        # Get next BUFR message
        ibufr = ec.codes_bufr_new_from_file(fptr)
        n = n + 1

        continue


"""
End read_bufr.py
""" 
