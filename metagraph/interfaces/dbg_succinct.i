%module metagraph
%include "std_string.i"
%include "std_vector.i"
%include "std_set.i"
%include "std_map.i"
%include "stdint.i"

%{
#include "dbg_succinct.hpp"
%}

%include "dbg_succinct.hpp"
