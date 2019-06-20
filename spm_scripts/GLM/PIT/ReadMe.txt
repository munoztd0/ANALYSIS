ReadMe PIT:


GLM-01: no modulators & durations = 0 (stick functions)
--> 4 basic contrasts (CSp-CSm, CSp-Base,  CSp-CSm&Base, grips)

GLM-02: no modulators & durations = 1
--> 4 basic contrasts (CSp-CSm, CSp-Base,  CSp-CSm&Base, grips)

GLM-03: 1st level modulators &  durations = 1 & Z-scored for fsl.txt
-->  9 contrasts (4 basics +CSm-Base + CSp*eff-CSm*eff + CSp*eff-Base*eff + CSp*eff-CSm*eff&Base*eff +CSm*eff-Base*eff)

GLM-04: 2nd level covariate & durations = 1
-->  9 contrasts (5 basic + 4*eff)

GLM-05: 2nd level covariate (diff of the Z-scored fsl.txt) & durations = 1
-->  9 contrasts (5 basic + 4*eff)