ReadMe PIT:


GLM-01: no modulators & durations = 0 (stick functions)
--> 5 basic contrasts (CSp-CSm, CSp-Base,  CSp-CSm&Base, grips, CSm-Base)

GLM-02: no modulators & durations = 1
--> 5 basic contrasts (CSp-CSm, CSp-Base,  CSp-CSm&Base, grips, CSm-Base)

GLM-03: 1st level modulators &  durations = 1 & Z-scored for fsl.txt
-->  9 contrasts (5 basic +CSm-Base + CSp*eff-CSm*eff + CSp*eff-Base*eff + CSp*eff-CSm*eff&Base*eff +CSm*eff-Base*eff)

GLM-04: 2nd level covariate (diff of the Z-scored fsl.txt) & durations = 1
-->  9 contrasts (5 basic + 4*eff)