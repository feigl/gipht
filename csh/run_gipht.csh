#!/bin/csh

matlab -nodisplay <<!
%giphtpath
run(strcat(getenv('HOME'),filesep,'gipht',filesep,'src',filesep','giphtpath'));
gipht
!


