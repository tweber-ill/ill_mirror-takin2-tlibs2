/**
 * tlibs2 -- swig interface
 * @author Tobias Weber
 * @date 4-jun-2020
 * @license see 'LICENSE' file
 */

%module tl2
%{
	#include "instr.h"
%}


%include "instr.h"

%template(FileInstrBaseD) tl2::FileInstrBase<double>;
