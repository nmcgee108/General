
import math

# We know that python uses IEEE double precision by default

i = prev = 3/2

     # Store the previous fl-pt here
     # (Because the "overlap" is caused by a loss
     # of precision, the values with the same 
     # reciprocal should be consecutive)

while i < 2:  

    i = math.nextafter(i, 2.0)   # Iterates to the next fl-pt number
    # i += 2**(-52)              # A hard-coded alternative for 3/2<i<2
    
          
    if 1/i == 1/prev:                       # If the recips. are equal:
                                            
        print(i.hex(), prev.hex())          # Print the hex representation
        print(f"{i:.20f}", f"{prev:.20f}")  # Print the normal representation
        
        break                               # Stop looking once we find a pair
        
    prev = i           # Store the new "previous" value

    