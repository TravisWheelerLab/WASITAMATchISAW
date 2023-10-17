#!/bin/bash
for m in mat/*; do o=${m%.mat}; ./lambda freqs.txt $m > out/${o#mat/}.out; done;
