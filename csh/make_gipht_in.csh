#!/bin/csh 

#foreach ING (*_gipht.in)
#foreach TAG (T344c_jul14best T115c_jul14best T72_jul14best T222_jul14best VT115_jul14best VT344_jul14best T451jul14best)

#set TAG = `echo $ING | sed 's/_gipht.in//'`
#hengill.geology.wisc.edu% ls -1 ../interferogram_lists/E*Biggs.lst | awk '{printf("#%s\n",substr($1,24,8))}'
foreach TAG (ERS_T072 ERS_T115 ERS_T179 ERS_T344 ERS_T451 ENV_T115 ENV_T344)

echo -n $TAG

cat MASTER5.ING | grep -iv track | grep -iv test >! tmp1.ing

cat tmp1.ing | grep -v ilist   >! tmp2.ing;     cat MASTER5.ING | grep ilist  | grep $TAG | sed 's/%//'    >> tmp2.ing; \mv tmp2.ing tmp1.ing
cat tmp1.ing | grep -v unitv_  >! tmp2.ing;     cat MASTER5.ING | grep unitv_ | grep $TAG | sed 's/%//'    >> tmp2.ing; \mv tmp2.ing tmp1.ing

cat tmp1.ing | grep -v nprocessors >! tmp2.ing ; echo 'nprocessors = 12'   >> tmp2.ing ; \mv tmp2.ing tmp1.ing
cat tmp1.ing | grep -v fitfun      >! tmp2.ing ; echo 'fitfun = funfit28'  >> tmp2.ing ; \mv tmp2.ing tmp1.ing

\mv -f tmp1.ing $TAG.ing

if (! -z $TAG.ing) then
    echo " created $TAG.ing"
else
    echo ' failed'
endif

end # loop over files

