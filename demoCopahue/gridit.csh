#!/bin/csh

matlab -nodisplay << EOF
load_quad
EOF

cat Quad_ASC_CSK_2013.75_2014_04_09_lonlatvlos.csv | sed 's/,/ /g' | awk '{print $0}' >! tmp.txt

blockmean -I0.004444/0.004445 `gmtinfo -I0.004444/0.004445 tmp.txt` tmp.txt \
| xyz2grd -I0.004444/0.004445 `gmtinfo -I0.004444/0.004445 tmp.txt` -GQuad_ASC_CSK_2013.75_2014_04_09_lonlatvlos.grd

matlab -nodisplay << EOF
load_quad
EOF

