


GLM-03: parametric modulators only on ANTs no classifier for denoise
GLM-04: parametric modulators only on ANTS + classifier for denies

GLM-05: parametric modulators on with durations rather than stick functions on ANTS + classifier denoise --> only liking as parametric modulator

GLM-06: parametric modulators on with durations rather than stick functions on ANTS no classifier denoise --> only liking as parametric modulator

GLM-07: parametric modulator 1; -1 for control vs reward with durations and on ANTs + classifier for denoise


GLM-08: parametric modulators with durations control vs reward with more accurate model of events of non-interest

GLM-09: parametric modulators with 2 second durations control vs reward with more accurate model of events of non-interest (to remove possible movement)

GLM-10: parametric modulators with durations control vs reward with more accurate model of events of non-interest ; we entered the estimated onset of the swall signal to remove movement

TODO/TOTRY

Using GLM-09 as guide
GLM-11: but on the swallowing cue; with duration (contrast taste vs tasteless --> do we get posterior insula?)
GLM-12: add paramentric modulators (attention add "SPM.Sess(ses).U(c).P(nc).orth  = 0;" (Ask david for the code to insert the modulator only if is non 0 and varies (e.g., the participant did not put all the time 50): onset liquid (moduluated by intensity,familiarity and liking)

Ev se disocccupata: glm-13 CON COVARIATA SECOND LIVELLO -> chiedi a David di mostrarti quale script da REWOD


