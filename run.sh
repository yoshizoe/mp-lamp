#!/bin/sh

#mpiexec -np  3 ./cont-lamp --item ./samples/cont_data/synth_100_10.data --pos  ./samples/cont_data/synth_100_10.class  --show_progress true --log true
mpiexec -np  3 ./cont-lamp --item ./samples/cont_data/synth_400_30.data --pos  ./samples/cont_data/synth_400_30.class  --show_progress true --log true
#mpiexec -np  3 ./cont-lamp --item ./samples/cont_data/synth_1000_20.data --pos  ./samples/cont_data/synth_1000_20.class  --show_progress true --log true

