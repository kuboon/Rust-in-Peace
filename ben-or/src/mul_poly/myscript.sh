#!/bin/bash
echo “Make c-progam 'time_measurement'“

gcc time_measurement.c -o time_measurement -lm  # ohne Optimierung
# gcc time_measurement.c -o time_measurement -lm -O2   # mit Optimierung

