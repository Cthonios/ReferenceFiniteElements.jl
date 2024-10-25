# Storage Types
We can setup ```ReferenceFE``` with several different backend for efficiency reasons.

For ```StaticArrays``` that are large than ~50 the compiler really starts to struggle,
so we can fall back to traditional arrays here. The only issue there is on GPUs. 
The data structures can not be transferred to GPUs for very high order elements.
