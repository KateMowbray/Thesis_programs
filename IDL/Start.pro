COMMON wkdirs, wkdir_global, tt, FRon, len
@ IDL/StartCFD.pro

; km: I'v e adjusted this slightly to also include the CMT specific routines I've written/modified

;.r IDL/getdata.pro  	    ; reads in data
;.r IDL/getndata.pro  	    ; reads in nonrelativistic data
.r getrdata.pro  	    ; reads in relativistic data
.r iniorbplot.pro	    ; plot initial array of positions
.r legend.pro   	    ; creates quick legend
.r mk_vector.pro	    ; makes pointy 3d arrows
;.r ODEINT.pro   	    ; ODE integrator
.r CMTODEINT.pro
.r orb__define.pro	    ; used to define 3d orbs
.r particletrack.pro    ; for 1 particle highlights track in 3d
.r CMTparticletrack.pro
.r CMTparticletrackcp.pro
.r CMTparticletrackfp.pro
.r CMTframegen
.r CMTerror
;.r IDL/pvels.pro    	    ; creates 3d vector plot
;.r IDL/remove.pro   	    ; removes elements from array
;.r IDL/regionzoom.pro	    ; like particletrack but zoomable
.r symbol_obj.pro	    ; ?
.r quickbifur.pro	    ; a program to (coarsely) look for zeros in 1d
.r StartLARE.pro
;.r ODEINT
.r fieldwrapping

Q0 = 1.60217646d-19 ; proton charge [C]
M0 = 9.10938188d-31 ; electron mass [kg]
kb = 1.3806503d-23  ; Boltzmann's constant [J/K]
wkdir_global="Data"

.r trilinear_3D tracefield_3D 
.r destaggerB
.r xplot3dJT
;.r FLINE.pro
;.r CMTFLINE.pro



