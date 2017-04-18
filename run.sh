#!/bin/sh

mpiexec -np 3 ./bin-lamp --item samples/sample_data/sample_item.csv --pos samples/sample_data/sample_expression_over1.csv --show_progress true --log true
