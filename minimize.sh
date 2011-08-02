#!/bin/bash

# A tool to reduce to only one conformation per ring
# In future versions of Frog, this might be unwanted.

export nl=`grep -n '$$$' $1 | head -n 1 | cut -d':' -f1`
echo $nl
export cmd="head -n $nl '$1' > '$1'.one && mv '$1'.one '$1'"
echo $cmd
eval $cmd
